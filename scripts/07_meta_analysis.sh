#!/usr/bin/env bash
# =============================================================================
# scripts/07_meta_analysis.sh  —  Step VII: MANTA meta-analysis across cohorts
# =============================================================================
# Combines the per-cohort summary files written by Step IV
# (<cohort_dir>/meta/<stratum>/<label>.mvsumstats.tsv.gz and <label>.rescor.tsv)
# for every cohort run listed in META_INPUTS, and writes one multivariate
# p-value per variant (P_manta_meta) to
#   OUTPUT_DIR/meta_analysis/<stratum>/meta_manta.tsv.gz
# Only summary files are read: no individual-level data leaves a cohort run.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

META="${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
if [ -f "${META}" ]; then
    while IFS='=' read -r key val; do
        [[ "${key}" =~ ^# ]] && continue; [[ -z "${key}" ]] && continue
        export "${key}=${val}"
    done < "${META}"
fi

MIN_COH="${META_MIN_COHORTS:-2}"
[ -n "${META_INPUTS:-}" ] || die "Step VII: META_INPUTS is empty. List the cohort OUTPUT_DIRs, comma-separated."
IFS=',' read -ra COHORT_DIRS <<< "${META_INPUTS}"
[ "${#COHORT_DIRS[@]}" -ge 2 ] || die "Step VII: META_INPUTS lists ${#COHORT_DIRS[@]} cohort(s); at least 2 are needed."

# ── 7A. Check every cohort run: folder, Step IV export, genome build ─────────
log_section "[Step VII] Checking cohort runs"
BUILDS=()
for d in "${COHORT_DIRS[@]}"; do
    [ -d "${d}" ] || die "Step VII: cohort folder not found: ${d}"
    [ -d "${d}/meta" ] || die "Step VII: no meta/ folder in ${d}. Run that cohort's Steps III-IV with PER_TRAIT_STATS=true first."
    md=$(ls "${d}"/*_run_metadata.txt 2>/dev/null | head -1 || true)
    b=""
    [ -n "${md}" ] && b=$(grep "^GENOME_BUILD=" "${md}" | tail -1 | cut -d= -f2-)
    BUILDS+=("${b:-unknown}")
    log_info "  ${d}  (genome build: ${b:-unknown})"
done
if [ "$(printf '%s\n' "${BUILDS[@]}" | grep -v '^unknown$' | sort -u | wc -l)" -gt 1 ]; then
    die "Step VII: cohorts use different genome builds (${BUILDS[*]}). Variants are matched by CHR:POS:A1:A2, so all cohorts must be on the same build."
fi

# ── 7B. One meta-analysis per stratum present in every cohort ────────────────
read -ra STRATA_ARR <<< "${STRATA:-combined}"
N_DONE=0
for STRATUM in "${STRATA_ARR[@]}"; do
    log_section "[Step VII] Stratum: ${STRATUM}"
    SS=(); RC=(); OK=true
    for d in "${COHORT_DIRS[@]}"; do
        s=( "${d}/meta/${STRATUM}"/*.mvsumstats.tsv.gz )
        r=( "${d}/meta/${STRATUM}"/*.rescor.tsv )
        if [ ! -f "${s[0]}" ] || [ ! -f "${r[0]}" ]; then
            log_warn "  ${d}: no Step IV export for stratum '${STRATUM}' — stratum skipped"
            OK=false; break
        fi
        if [ "${#s[@]}" -ne 1 ] || [ "${#r[@]}" -ne 1 ]; then
            die "Step VII: expected one *.mvsumstats.tsv.gz and one *.rescor.tsv in ${d}/meta/${STRATUM}"
        fi
        SS+=("${s[0]}"); RC+=("${r[0]}")
    done
    [ "${OK}" = "true" ] || continue

    OUT_DIR="${OUTPUT_DIR}/meta_analysis/${STRATUM}"
    mkdir -p "${OUT_DIR}"
    {
        echo -e "cohort_dir\tsumstats\trescor"
        for i in "${!COHORT_DIRS[@]}"; do
            echo -e "${COHORT_DIRS[$i]}\t${SS[$i]}\t${RC[$i]}"
        done
    } > "${OUT_DIR}/meta_inputs.tsv"

    mkdir -p "${OUTPUT_DIR}/logs"
    run_cmd "  meta_manta.R on ${#SS[@]} cohorts [${STRATUM}]" \
        Rscript "${SCRIPT_DIR}/scripts/meta_manta.R" \
        --sumstats    "$(IFS=,; echo "${SS[*]}")" \
        --rescor      "$(IFS=,; echo "${RC[*]}")" \
        --min-cohorts "${MIN_COH}" \
        --out         "${OUT_DIR}/meta_manta.tsv.gz" \
        2>&1 | tee -a "${OUTPUT_DIR}/logs/meta_analysis_${STRATUM}.log"

    STRATUM_UC=$(echo "${STRATUM}" | tr 'a-z' 'A-Z')
    echo "META_RESULT_${STRATUM_UC}=${OUT_DIR}/meta_manta.tsv.gz" >> "${META}"
    log_ok "  [${STRATUM}] ${OUT_DIR}/meta_manta.tsv.gz"
    N_DONE=$((N_DONE + 1))
done

[ "${N_DONE}" -gt 0 ] || die "Step VII: no stratum could be meta-analysed (see warnings above)."
record_stat "meta_strata" "${N_DONE}"
