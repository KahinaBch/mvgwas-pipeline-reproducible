#!/usr/bin/env bash
# =============================================================================
# scripts/03_run_pipeline.sh  —  Step III: mvgwas-nf execution + merge + QC
# =============================================================================
# Splits the VCF by chromosome, runs mvgwas-nf per chromosome (locally in
# parallel or via SLURM), merges the results, and runs coherence checks.
#
# Reads extra run-metadata written by Step II:
#   PHENOTYPE_FILTERED, COVARIATE_FILTERED, SAMPLES_INTERSECTION, VCF_CHR_PREFIX
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

# ── Re-read metadata from Step II ────────────────────────────────────────────
META="${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
if [ -f "${META}" ]; then
    while IFS='=' read -r key val; do
        [[ "${key}" =~ ^# ]] && continue
        [[ -z "${key}" ]] && continue
        export "${key}=${val}"
    done < "${META}"
fi

PHENO_FINAL="${PHENOTYPE_FILTERED:-${OUTPUT_DIR}/inputs/phenotype_filtered.tsv}"
COV_FINAL="${COVARIATE_FILTERED:-${OUTPUT_DIR}/inputs/covariate_filtered.tsv}"
SAMPLES_FILE="${SAMPLES_INTERSECTION:-${OUTPUT_DIR}/inputs/samples_intersection.txt}"

require_file "${PHENO_FINAL}"    "Filtered phenotype"
require_file "${COV_FINAL}"      "Filtered covariate"
require_file "${SAMPLES_FILE}"   "Sample intersection"

N_SAMPLES=$(wc -l < "${SAMPLES_FILE}")
log_info "  Analysis sample size: ${N_SAMPLES}"

# Directories
CHR_VCF_DIR="${OUTPUT_DIR}/chr_vcfs"
CHR_RESULTS_DIR="${OUTPUT_DIR}/chr_results"
LOGS_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${CHR_VCF_DIR}" "${CHR_RESULTS_DIR}" "${LOGS_DIR}"

# Container flag (exported by Step I, or fall back from config)
NXF_CONTAINER_FLAG="${NXF_CONTAINER_FLAG:-}"
case "${CONTAINER:-singularity}" in
    singularity) NXF_CONTAINER_FLAG="${NXF_CONTAINER_FLAG:--with-singularity}" ;;
    docker)      NXF_CONTAINER_FLAG="${NXF_CONTAINER_FLAG:--with-docker}"      ;;
    none)        NXF_CONTAINER_FLAG="" ;;
esac

# Nextflow version prefix
NXF_PREFIX=""
[ -n "${NXF_VER:-}" ] && NXF_PREFIX="NXF_VER=${NXF_VER} "

# Custom nextflow config
NXF_CONFIG_ARG=""
if [ -n "${NXF_CONFIG:-}" ] && [ -f "${NXF_CONFIG}" ]; then
    NXF_CONFIG_ARG="-c ${NXF_CONFIG}"
elif [ -f "${SCRIPT_DIR}/nextflow/nextflow.config" ]; then
    NXF_CONFIG_ARG="-c ${SCRIPT_DIR}/nextflow/nextflow.config"
fi

# =============================================================================
# 3A. Split VCF by chromosome (with sample subsetting)
# =============================================================================
log_section "[Step III] Splitting VCF by chromosome"

read -ra CHR_ARRAY <<< "${CHROMOSOMES:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"

for CHR in "${CHR_ARRAY[@]}"; do
    CHR_VCF="${CHR_VCF_DIR}/chr${CHR}.vcf.gz"
    # Determine the chromosome label to query from the VCF
    CHR_LABEL="${VCF_CHR_PREFIX:-}${CHR}"

    if [ -f "${CHR_VCF}" ] && [ -f "${CHR_VCF}.tbi" ]; then
        log_debug "  chr${CHR}: VCF already exists, skipping extraction"
        continue
    fi

    log_info "  Extracting chr${CHR} from ${VCF_FILE}..."
    if [ "${DRY_RUN}" = "false" ]; then
        bcftools view \
            -r "${CHR_LABEL}" \
            -S "${SAMPLES_FILE}" \
            --force-samples \
            "${VCF_FILE}" \
            -Oz -o "${CHR_VCF}" \
            2>> "${LOGS_DIR}/bcftools_split.log"
        bcftools index -t "${CHR_VCF}"
        N_VARS=$(bcftools index -n "${CHR_VCF}" 2>/dev/null || echo "?")
        log_ok "    chr${CHR}: ${N_VARS} variants → ${CHR_VCF}"
    else
        log_info "  [DRY-RUN] Would extract chr${CHR}"
    fi
done

log_ok "  VCF splitting complete."

# =============================================================================
# 3B. Run mvgwas-nf per chromosome
# =============================================================================
log_section "[Step III] Running mvgwas-nf per chromosome"

EXEC_MODE=$(detect_environment)
log_info "  Execution mode: ${EXEC_MODE}  (USE_SLURM=${USE_SLURM:-false})"
log_info "  Chromosomes   : ${CHR_ARRAY[*]}"

# ── Helper: build nextflow command for one chromosome ────────────────────────
build_nxf_cmd() {
    local chr="$1"
    local chr_vcf="${CHR_VCF_DIR}/chr${chr}.vcf.gz"
    local chr_out="${CHR_RESULTS_DIR}/chr${chr}"
    mkdir -p "${chr_out}"
    echo "${NXF_PREFIX}nextflow run ${PIPELINE_DIR}/mvgwas.nf" \
        "--geno ${chr_vcf}" \
        "--pheno ${PHENO_FINAL}" \
        "--cov ${COV_FINAL}" \
        "--dir ${chr_out}" \
        "--out mvgwas_chr${chr}.tsv" \
        "--l ${CHUNK_SIZE:-500}" \
        "--t ${TRANSFORMATION:-none}" \
        "--i ${INTERACTION:-none}" \
        "--ng ${MIN_INDIVIDUALS:-10}" \
        "${NXF_CONTAINER_FLAG}" \
        "${NXF_CONFIG_ARG}" \
        "-work-dir ${chr_out}/.nf_work" \
        "-resume"
}

if [ "${EXEC_MODE}" = "slurm" ]; then
    # ── SLURM job array approach ─────────────────────────────────────────────
    log_info "  Submitting SLURM job array..."

    # Build array spec from CHROMOSOMES (e.g. "1 2 3" → "1,2,3")
    ARRAY_SPEC=$(echo "${CHROMOSOMES:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}" | tr ' ' ',')

    # Write the per-chromosome worker script
    CHR_WORKER="${LOGS_DIR}/slurm_chr_worker.sh"
    cat > "${CHR_WORKER}" <<EOSLURM
#!/bin/bash
#SBATCH --job-name=${RUN_NAME}_mvgwas
#SBATCH --output=${LOGS_DIR}/chr_%a.out
#SBATCH --error=${LOGS_DIR}/chr_%a.err
#SBATCH --time=${SLURM_TIME:-24:00:00}
#SBATCH --mem=${SLURM_MEM:-32G}
#SBATCH --cpus-per-task=${SLURM_CPUS:-4}
#SBATCH --array=${ARRAY_SPEC}
#SBATCH --partition=${SLURM_PARTITION:-genoa}
$([ -n "${SLURM_MAIL:-}" ] && echo "#SBATCH --mail-type=FAIL" && echo "#SBATCH --mail-user=${SLURM_MAIL:-}" || echo "")
${SLURM_EXTRA_ARGS:-}

CHR=\${SLURM_ARRAY_TASK_ID}
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/utils.sh"
export PIPELINE_LOG="${LOGS_DIR}/chr_\${CHR}.log"
export OUTPUT_DIR="${OUTPUT_DIR}"
export DRY_RUN="${DRY_RUN}"
export VERBOSE="${VERBOSE}"
export RUN_NAME="${RUN_NAME}"

echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Starting chr\${CHR}"
$(build_nxf_cmd "\${CHR}") 2>&1 | tee -a "\${PIPELINE_LOG}"
echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Finished chr\${CHR}"
EOSLURM

    ARRAY_JOB_ID=$(sbatch --parsable "${CHR_WORKER}")
    log_ok "  Submitted SLURM array job: ${ARRAY_JOB_ID} (22 tasks)"

    # Submit merge job with dependency
    MERGE_JOB_ID=$(sbatch --parsable \
        --dependency=afterok:${ARRAY_JOB_ID} \
        --job-name="${RUN_NAME}_merge" \
        --output="${LOGS_DIR}/merge.out" \
        --error="${LOGS_DIR}/merge.err" \
        --partition="${SLURM_PARTITION:-genoa}" \
        --time="${SLURM_MERGE_TIME:-02:00:00}" \
        --mem="${SLURM_MERGE_MEM:-16G}" \
        --cpus-per-task="${SLURM_MERGE_CPUS:-2}" \
        --wrap="bash ${SCRIPT_DIR}/scripts/03_merge_results.sh")
    log_ok "  Submitted merge job: ${MERGE_JOB_ID} (depends on ${ARRAY_JOB_ID})"
    log_info "  Monitor: squeue -j ${ARRAY_JOB_ID},${MERGE_JOB_ID}"
    log_warn "  SLURM mode: run.sh will now exit. Re-run with --steps 4,5 --resume after jobs complete."
    record_stat "slurm_array_job_id" "${ARRAY_JOB_ID}"
    record_stat "slurm_merge_job_id" "${MERGE_JOB_ID}"
    exit 0

else
    # ── Local parallel execution ─────────────────────────────────────────────
    LOCAL_PARALLEL="${LOCAL_PARALLEL:-4}"
    log_info "  Running up to ${LOCAL_PARALLEL} chromosomes in parallel"

    BG_PIDS=()
    ACTIVE=0

    for CHR in "${CHR_ARRAY[@]}"; do
        CHR_OUT="${CHR_RESULTS_DIR}/chr${CHR}"
        RESULT_FILE="${CHR_OUT}/result/mvgwas_chr${CHR}.tsv"

        if [ -f "${RESULT_FILE}" ] && [ -s "${RESULT_FILE}" ]; then
            log_info "  chr${CHR}: result already exists ($(wc -l < "${RESULT_FILE}") lines) — skipping"
            continue
        fi

        NXF_CMD=$(build_nxf_cmd "${CHR}")
        log_info "  Launching chr${CHR}..."
        log_debug "    CMD: ${NXF_CMD}"

        if [ "${DRY_RUN}" = "false" ]; then
            (
                mkdir -p "${CHR_OUT}"
                cd "${CHR_OUT}"
                eval "${NXF_CMD}" \
                    >> "${LOGS_DIR}/chr${CHR}.log" 2>&1
                echo "EXIT:$?" >> "${LOGS_DIR}/chr${CHR}.log"
            ) &
            BG_PIDS+=($!)
            ACTIVE=$(( ACTIVE + 1 ))

            # Throttle parallelism
            if [ "${LOCAL_PARALLEL}" -gt 0 ] && [ "${ACTIVE}" -ge "${LOCAL_PARALLEL}" ]; then
                wait "${BG_PIDS[0]}" || log_warn "  A chromosome job may have failed; check logs"
                BG_PIDS=("${BG_PIDS[@]:1}")
                ACTIVE=$(( ACTIVE - 1 ))
            fi
        else
            log_info "  [DRY-RUN] Would run: ${NXF_CMD}"
        fi
    done

    # Wait for remaining jobs
    if [ "${DRY_RUN}" = "false" ] && [ ${#BG_PIDS[@]} -gt 0 ]; then
        log_info "  Waiting for final chromosome jobs..."
        FAILED=0
        for pid in "${BG_PIDS[@]}"; do
            wait "${pid}" || (( FAILED++ )) || true
        done
        if [ "${FAILED}" -gt 0 ]; then
            die "${FAILED} chromosome jobs failed. Check logs in ${LOGS_DIR}/"
        fi
        log_ok "  All chromosome jobs finished."
    fi
fi

# =============================================================================
# 3C. Merge chromosome results
# =============================================================================
bash "${SCRIPT_DIR}/scripts/03_merge_results.sh"
