#!/usr/bin/env bash
# =============================================================================
# run.sh  —  mvgwas-pipeline: Automated multivariate GWAS runner
# =============================================================================
#
# Six-step pipeline:
#   I.   Environment setup and pipeline self-test (skippable)
#   II.  Input format validation and individual harmonisation
#   III. mvgwas-nf execution (per-chromosome parallel) + output QC
#   IV.  Top-SNP extraction, rsID annotation, Manhattan/QQ/regional plots
#   V.   Enrichment analysis: GWAS Catalog mapping + MAGMA gene/gene-set test
#          (runs by default; skip with --skip-enrichment or SKIP_ENRICHMENT=true)
#   VI.  Comprehensive final report
#
# Usage:
#   bash run.sh --config <config.conf> [options]
#
# Options:
#   --config  <file>            Configuration file (default: config.conf)
#   --skip-test                 Skip Step I (environment + pipeline test)
#   --build   GRCh37|GRCh38    Override genome build for rsID annotation
#   --steps   <list>            Run only these steps, e.g. "2,3,4" (default: 1-6)
#   --resume                    Skip steps that already have a checkpoint
#   --dry-run                   Print commands without executing
#   --verbose                   Enable debug-level logging
#   --skip-enrichment           Skip Step V (GWAS Catalog + MAGMA)
#   --sex-stratified            Run sex-stratified analysis (male/female/combined)
#   --help                      Show this message and exit
#
# Examples:
#   bash run.sh --config my_run.conf
#   bash run.sh --config my_run.conf --skip-test --resume
#   bash run.sh --config my_run.conf --steps 4,5 --build GRCh38
# =============================================================================
set -euo pipefail

PIPELINE_START_TS=$(date +%s)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Source libraries ──────────────────────────────────────────────────────────
# shellcheck source=lib/logger.sh
source "${SCRIPT_DIR}/lib/logger.sh"
# shellcheck source=lib/utils.sh
source "${SCRIPT_DIR}/lib/utils.sh"

# ── Default option values ─────────────────────────────────────────────────────
CONFIG_FILE="${SCRIPT_DIR}/config.conf"
SKIP_TEST=false
BUILD_OVERRIDE=""
STEPS_TO_RUN="1,2,3,4,5,6"
RESUME=false
DRY_RUN=false
VERBOSE=false
SKIP_ENRICHMENT=false

# ── Parse command-line arguments ──────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)    CONFIG_FILE="$2";      shift 2 ;;
        --skip-test) SKIP_TEST=true;        shift   ;;
        --build)     BUILD_OVERRIDE="$2";   shift 2 ;;
        --steps)     STEPS_TO_RUN="$2";     shift 2 ;;
        --resume)    RESUME=true;           shift   ;;
        --dry-run)   DRY_RUN=true;          shift   ;;
        --verbose)           VERBOSE=true;          shift   ;;
        --sex-stratified)    SEX_STRATIFIED=true;   shift   ;;
        --skip-enrichment)   SKIP_ENRICHMENT=true;  shift   ;;
        --no-enrichment)     SKIP_ENRICHMENT=true;  shift   ;;
        --help|-h)
            grep '^#' "${BASH_SOURCE[0]}" | head -40 | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Load configuration file ───────────────────────────────────────────────────
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Configuration file not found: ${CONFIG_FILE}"
    echo "Copy config.conf.template to config.conf and fill in the paths."
    exit 1
fi
# shellcheck source=config.conf.template
source "${CONFIG_FILE}"

# Apply overrides
[ -n "${BUILD_OVERRIDE}" ] && GENOME_BUILD="${BUILD_OVERRIDE}"
export DRY_RUN VERBOSE SEX_STRATIFIED SKIP_ENRICHMENT

# ── Validate mandatory config keys ───────────────────────────────────────────
_check_config() {
    local key="$1"; local val="${!key:-}"
    if [ -z "${val}" ]; then
        echo "ERROR: Required config key '${key}' is not set in ${CONFIG_FILE}."
        exit 1
    fi
}
_check_config VCF_FILE
_check_config PHENOTYPE_FILE
_check_config COVARIATE_FILE
_check_config OUTPUT_DIR
_check_config PIPELINE_DIR
_check_config RUN_NAME

# ── Set up output directory and log file ─────────────────────────────────────
mkdir -p "${OUTPUT_DIR}"
PIPELINE_LOG="${OUTPUT_DIR}/${RUN_NAME}_pipeline.log"
export PIPELINE_LOG OUTPUT_DIR

# ── Banner ────────────────────────────────────────────────────────────────────
log_banner "mvgwas-pipeline  |  Run: ${RUN_NAME}  |  $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Config     : ${CONFIG_FILE}"
log_info "Output     : ${OUTPUT_DIR}"
log_info "Log        : ${PIPELINE_LOG}"
log_info "Genome     : ${GENOME_BUILD}"
log_info "Steps      : ${STEPS_TO_RUN}"
log_info "Resume     : ${RESUME}"
log_info "Dry-run    : ${DRY_RUN}"
log_info "SLURM      : ${USE_SLURM:-false}"
log_info ""

# Write run metadata
{
    echo "# mvgwas-pipeline run metadata"
    echo "RUN_NAME=${RUN_NAME}"
    echo "START_TIME=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "CONFIG_FILE=${CONFIG_FILE}"
    echo "GENOME_BUILD=${GENOME_BUILD}"
    echo "VCF_FILE=${VCF_FILE}"
    echo "PHENOTYPE_FILE=${PHENOTYPE_FILE}"
    echo "COVARIATE_FILE=${COVARIATE_FILE}"
    echo "PIPELINE_DIR=${PIPELINE_DIR}"
    echo "SEX_STRATIFIED=${SEX_STRATIFIED:-false}"
    echo "SEX_COL=${SEX_COL:-}"
    echo "SKIP_ENRICHMENT=${SKIP_ENRICHMENT:-false}"
} > "${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"

# Initialise run_stats.tsv
echo -e "key\tvalue" > "${OUTPUT_DIR}/run_stats.tsv"
record_stat "run_name"    "${RUN_NAME}"
record_stat "start_time"  "$(date '+%Y-%m-%d %H:%M:%S')"
record_stat "genome_build" "${GENOME_BUILD}"

# ── Helper: should_run_step ───────────────────────────────────────────────────
should_run_step() {
    local n="$1"
    if [[ ",${STEPS_TO_RUN}," != *",${n},"* ]]; then
        log_debug "Step ${n} not in --steps list; skipping."
        return 1
    fi
    if [ "${RESUME}" = "true" ] && checkpoint_exists "${n}"; then
        local ckpt_ts; ckpt_ts=$(cat "${OUTPUT_DIR}/.checkpoints/step_${n}.done")
        log_info "Step ${n}: checkpoint found (${ckpt_ts}) — skipping (--resume)."
        return 1
    fi
    return 0
}

# ── STEP I: Environment setup and pipeline test ───────────────────────────────
if should_run_step 1; then
    if [ "${SKIP_TEST}" = "true" ]; then
        log_step 1 "Environment setup (TEST SKIPPED — --skip-test flag set)"
        bash "${SCRIPT_DIR}/scripts/01_setup_test.sh" --no-test
    else
        log_step 1 "Environment setup and pipeline self-test"
        bash "${SCRIPT_DIR}/scripts/01_setup_test.sh"
    fi
    checkpoint_set 1
    log_milestone "Step I complete — environment verified"
    record_stat "step1_completed" "$(date '+%Y-%m-%d %H:%M:%S')"
fi

# ── STEP II: Input validation and individual harmonisation ────────────────────
if should_run_step 2; then
    log_step 2 "Input format validation and individual harmonisation"
    bash "${SCRIPT_DIR}/scripts/02_check_inputs.sh"
    checkpoint_set 2
    log_milestone "Step II complete — inputs validated and harmonised"
    record_stat "step2_completed" "$(date '+%Y-%m-%d %H:%M:%S')"
fi

# Read STRATA from metadata written by Step II (default: combined)
META_FILE="${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
if [ -f "${META_FILE}" ]; then
    _STRATA_VAL=$(grep "^STRATA=" "${META_FILE}" | tail -1 | cut -d= -f2-)
    [ -n "${_STRATA_VAL}" ] && STRATA="${_STRATA_VAL}"
    # Export all metadata keys so subscripts can resolve stratum-specific paths
    while IFS='=' read -r _K _V; do
        [[ "${_K}" =~ ^[A-Z_][A-Z0-9_]*$ ]] && export "${_K}=${_V}"
    done < "${META_FILE}"
fi
STRATA="${STRATA:-combined}"
export STRATA
log_info "Active strata: ${STRATA}"

# ── STEP III: Pipeline execution and output QC ────────────────────────────────
if should_run_step 3; then
    log_step 3 "mvgwas-nf execution (chromosome-parallel) and output QC"
    bash "${SCRIPT_DIR}/scripts/03_run_pipeline.sh"
    # Merge is called from inside 03_run_pipeline.sh after all chr jobs finish
    checkpoint_set 3
    log_milestone "Step III complete — GWAS complete, results merged and QC'd"
    record_stat "step3_completed" "$(date '+%Y-%m-%d %H:%M:%S')"
fi

# ── STEP IV: Post-processing (top SNPs, annotation, plots) ───────────────────
if should_run_step 4; then
    log_step 4 "Top-SNP extraction, rsID annotation, and visualisation"
    bash "${SCRIPT_DIR}/scripts/04_postprocess.sh"
    checkpoint_set 4
    log_milestone "Step IV complete — top SNPs annotated, plots generated"
    record_stat "step4_completed" "$(date '+%Y-%m-%d %H:%M:%S')"
fi

# ── STEP V: Enrichment analysis (GWAS Catalog + MAGMA) ───────────────────────
if should_run_step 5; then
    # Allow config or CLI flag to skip enrichment
    # Config file is already sourced; SKIP_ENRICHMENT may have been set there.
    if [ "${SKIP_ENRICHMENT:-false}" = "true" ]; then
        log_step 5 "Enrichment analysis — SKIPPED (SKIP_ENRICHMENT=true)"
    else
        log_step 5 "Enrichment analysis: GWAS Catalog mapping + MAGMA"
        read -ra STRATA_ARR_V <<< "${STRATA:-combined}"
        mkdir -p "${OUTPUT_DIR}/logs"

        for STRATUM in "${STRATA_ARR_V[@]}"; do
            log_section "  [Step V] Stratum: ${STRATUM}"
            STRATUM_UC=$(echo "${STRATUM}" | tr 'a-z' 'A-Z')

            # ── 5A. GWAS Catalog mapping ──────────────────────────────────────
            # Resolve annotated (rsID) top-SNPs file for this stratum
            if [ "${STRATUM}" = "combined" ]; then
                eval "_TOP_FILE=\"\${ANNOTATED_SNPS:-\${TOP_SNPS:-}}\""
            else
                eval "_TOP_FILE=\"\${ANNOTATED_SNPS_${STRATUM_UC}:-\${TOP_SNPS_${STRATUM_UC}:-}}\""
            fi

            CATALOG_DIR="${OUTPUT_DIR}/gwas_catalog/${STRATUM}"
            mkdir -p "${CATALOG_DIR}"

            if [ -n "${_TOP_FILE:-}" ] && [ -f "${_TOP_FILE}" ]; then
                log_info "  [5A] GWAS Catalog mapping [${STRATUM}]"
                python3 "${SCRIPT_DIR}/scripts/05_gwas_catalog.py" \
                    --top-snps   "${_TOP_FILE}" \
                    --output-dir "${CATALOG_DIR}" \
                    --run-name   "${RUN_NAME}" \
                    --stratum    "${STRATUM}" \
                    --pval-col   "${PVAL_COL:-P_manta}" \
                    --top-n      "${GWAS_CATALOG_TOP_N:-20}" \
                    --clump-kb   "${GWAS_CATALOG_CLUMP_KB:-500}" \
                    --window-kb  "${GWAS_CATALOG_WINDOW_KB:-500}" \
                    2>&1 | tee -a "${OUTPUT_DIR}/logs/gwas_catalog_${STRATUM}.log"

                echo "GWAS_CATALOG_REPORT_${STRATUM_UC}=${CATALOG_DIR}/gwas_catalog_${STRATUM}.md" \
                    >> "${META_FILE}"
                log_ok "  [5A] GWAS Catalog done → ${CATALOG_DIR}/gwas_catalog_${STRATUM}.md"
            else
                log_warn "  [5A] No annotated SNPs file for stratum '${STRATUM}' — skipping GWAS Catalog"
                log_warn "       Run Step IV first, or check that rsID annotation succeeded."
            fi

            # ── 5B. MAGMA ─────────────────────────────────────────────────────
            log_info "  [5B] MAGMA gene/gene-set enrichment [${STRATUM}]"
            bash "${SCRIPT_DIR}/scripts/05_magma.sh" "${STRATUM}" \
                2>&1 | tee -a "${OUTPUT_DIR}/logs/magma_${STRATUM}.log"
        done
    fi
    checkpoint_set 5
    log_milestone "Step V complete — enrichment analyses done"
    record_stat "step5_completed" "$(date '+%Y-%m-%d %H:%M:%S')"
fi

# ── STEP VI: Report generation ────────────────────────────────────────────────
if should_run_step 6; then
    log_step 6 "Comprehensive report generation"
    python3 "${SCRIPT_DIR}/scripts/05_generate_report.py" \
        --output-dir    "${OUTPUT_DIR}" \
        --run-name      "${RUN_NAME}" \
        --genome-build  "${GENOME_BUILD}" \
        --pval-col      "${PVAL_COL:-P_manta}" \
        --gw-thresh     "${GW_THRESH:-5e-8}" \
        --sug-thresh    "${SUG_THRESH:-1e-5}" \
        --config-file   "${CONFIG_FILE}" \
        --log-file      "${PIPELINE_LOG}"
    checkpoint_set 6
    log_milestone "Step VI complete — report written"
    record_stat "step6_completed" "$(date '+%Y-%m-%d %H:%M:%S')"
fi

# ── Final summary ─────────────────────────────────────────────────────────────
record_stat "end_time"    "$(date '+%Y-%m-%d %H:%M:%S')"
record_stat "elapsed"     "$(elapsed)"

log_banner "Pipeline complete!  Run: ${RUN_NAME}  |  Elapsed: $(elapsed)"
log_info "Results    : ${OUTPUT_DIR}/"
log_info "Report     : ${OUTPUT_DIR}/${RUN_NAME}_report.md"
log_info "GWAS Cat.  : ${OUTPUT_DIR}/gwas_catalog/"
log_info "MAGMA      : ${OUTPUT_DIR}/magma/"
log_info "Log        : ${PIPELINE_LOG}"
log_info ""
log_ok "All done."
