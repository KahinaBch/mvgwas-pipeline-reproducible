#!/usr/bin/env bash
# =============================================================================
# scripts/03_merge_results.sh  —  Merge per-chromosome mvgwas results + QC
# =============================================================================
# This script is either called directly from 03_run_pipeline.sh (local mode)
# or submitted as a dependent SLURM job (SLURM mode). It can also be invoked
# standalone after a SLURM run completes.
# =============================================================================
set -euo pipefail

# If run as a standalone SLURM job we need to re-source config
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

# Re-read metadata
META="${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
if [ -f "${META}" ]; then
    while IFS='=' read -r key val; do
        [[ "${key}" =~ ^# ]] && continue; [[ -z "${key}" ]] && continue
        export "${key}=${val}"
    done < "${META}"
fi

log_section "[Step III] Merging chromosome results"

read -ra CHR_ARRAY <<< "${CHROMOSOMES:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
CHR_RESULTS_DIR="${OUTPUT_DIR}/chr_results"
RESULTS_DIR="${OUTPUT_DIR}/results"
QC_DIR="${OUTPUT_DIR}/qc"
LOGS_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${RESULTS_DIR}" "${QC_DIR}"

MERGED="${RESULTS_DIR}/mvgwas_merged.tsv"
QC_REPORT="${QC_DIR}/merge_qc_report.txt"

# ── Locate per-chromosome result files ───────────────────────────────────────
declare -A CHR_FILE
MISSING_CHRS=()
PRESENT_CHRS=()

for CHR in "${CHR_ARRAY[@]}"; do
    # Try several naming patterns produced by different nextflow versions
    for candidate in \
        "${CHR_RESULTS_DIR}/chr${CHR}/result/mvgwas_chr${CHR}.tsv" \
        "${CHR_RESULTS_DIR}/chr${CHR}/mvgwas_chr${CHR}.tsv" \
        "${CHR_RESULTS_DIR}/chr${CHR}/results/mvgwas_chr${CHR}.tsv"; do
        if [ -f "${candidate}" ] && [ -s "${candidate}" ]; then
            CHR_FILE[$CHR]="${candidate}"
            PRESENT_CHRS+=("${CHR}")
            break
        fi
    done
    if [ -z "${CHR_FILE[$CHR]+x}" ]; then
        MISSING_CHRS+=("${CHR}")
        log_warn "  chr${CHR}: NO result file found"
    fi
done

N_PRESENT=${#PRESENT_CHRS[@]}
N_MISSING=${#MISSING_CHRS[@]}

log_info "  Present: ${N_PRESENT} chromosomes"
[ "${N_MISSING}" -gt 0 ] && log_warn "  Missing: ${N_MISSING} (${MISSING_CHRS[*]})"

if [ "${N_PRESENT}" -eq 0 ]; then
    die "No chromosome results found. Pipeline may not have run yet."
fi

if [ "${DRY_RUN}" = "false" ]; then
    # ── Merge: header once, then data rows from all present chromosomes ───────
    log_info "  Building merged file: ${MERGED}"
    FIRST_CHR="${PRESENT_CHRS[0]}"
    head -1 "${CHR_FILE[$FIRST_CHR]}" > "${MERGED}"

    TOTAL_ROWS=0
    declare -A CHR_ROWS

    for CHR in "${PRESENT_CHRS[@]}"; do
        NROWS=$(awk 'NR>1' "${CHR_FILE[$CHR]}" | wc -l)
        awk 'NR>1' "${CHR_FILE[$CHR]}" >> "${MERGED}"
        CHR_ROWS[$CHR]=${NROWS}
        TOTAL_ROWS=$(( TOTAL_ROWS + NROWS ))
        log_info "    chr${CHR}: ${NROWS} associations"
    done

    log_ok "  Merged: $(fmt_num ${TOTAL_ROWS}) total associations → ${MERGED}"
    record_stat "total_associations" "${TOTAL_ROWS}"
    record_stat "chromosomes_present" "${N_PRESENT}"
    record_stat "chromosomes_missing" "${N_MISSING}"

    # ── QC: detect outlier chromosomes + column check ─────────────────────────
    log_section "[Step III] QC checks on merged results"
    NCOLS=$(awk -F'\t' 'NR==1{print NF}' "${MERGED}")
    N_HEADER_ROWS=1
    N_TOTAL_LINES=$(wc -l < "${MERGED}")
    log_info "  Columns          : ${NCOLS}"
    log_info "  Total lines      : ${N_TOTAL_LINES}  (${TOTAL_ROWS} data rows)"

    # Check that PVAL_COL exists in header
    PVAL_COL="${PVAL_COL:-P_manta}"
    HEADER=$(head -1 "${MERGED}")
    if echo "${HEADER}" | grep -qw "${PVAL_COL}"; then
        log_ok "  P-value column '${PVAL_COL}' found"
    else
        log_warn "  P-value column '${PVAL_COL}' NOT found in header: ${HEADER}"
        log_warn "  Available columns: $(head -1 "${MERGED}" | tr '\t' '\n' | head -20 | tr '\n' ' ')"
    fi

    # Check for required columns
    for COL in CHR POS ID REF ALT; do
        if echo "${HEADER}" | tr '\t' '\n' | grep -qi "^${COL}$"; then
            log_ok "  Column '${COL}' found"
        else
            log_warn "  Column '${COL}' not found (might use different name)"
        fi
    done

    # Average rows per chromosome (detect very small chrs)
    AVG=$(( TOTAL_ROWS / N_PRESENT ))
    log_info "  Mean variants per chromosome: ${AVG}"
    for CHR in "${PRESENT_CHRS[@]}"; do
        RATIO=$(( CHR_ROWS[$CHR] * 100 / (AVG + 1) ))
        if [ "${RATIO}" -lt 20 ]; then
            log_warn "  chr${CHR} has only ${CHR_ROWS[$CHR]} rows (${RATIO}% of mean) — potential issue"
        fi
    done

    # ── Write QC report ───────────────────────────────────────────────────────
    {
        echo "=== mvgwas-pipeline — Merge QC Report ==="
        echo "Run name    : ${RUN_NAME}"
        echo "Timestamp   : $(date '+%Y-%m-%d %H:%M:%S')"
        echo ""
        echo "--- Input ---"
        echo "Chromosomes requested : ${CHR_ARRAY[*]}"
        echo "Chromosomes present   : ${PRESENT_CHRS[*]}"
        echo "Chromosomes missing   : ${MISSING_CHRS[*]:-none}"
        echo ""
        echo "--- Output ---"
        echo "Merged file : ${MERGED}"
        echo "Total rows  : ${TOTAL_ROWS}"
        echo "Columns     : ${NCOLS}"
        echo "Header      : ${HEADER}"
        echo ""
        echo "--- Per-chromosome row counts ---"
        for CHR in "${CHR_ARRAY[@]}"; do
            if [ -n "${CHR_FILE[$CHR]+x}" ]; then
                echo "  chr${CHR}: ${CHR_ROWS[$CHR]}"
            else
                echo "  chr${CHR}: MISSING"
            fi
        done
        echo ""
        echo "--- Status ---"
        [ "${N_MISSING}" -eq 0 ] && echo "PASS: all chromosomes present" || echo "WARN: ${N_MISSING} chromosomes missing"
    } > "${QC_REPORT}"
    log_ok "  QC report: ${QC_REPORT}"

    # ── Export merged path to metadata ────────────────────────────────────────
    echo "MERGED_RESULTS=${MERGED}" >> "${META}"
    echo "MERGE_QC_REPORT=${QC_REPORT}" >> "${META}"
    echo "TOTAL_ASSOCIATIONS=${TOTAL_ROWS}" >> "${META}"
    echo "HEADER_COLUMNS=${HEADER}" >> "${META}"

    log_milestone "Step III complete — merged results available at ${MERGED}"
else
    log_info "  [DRY-RUN] Would merge ${N_PRESENT} chromosome files → ${MERGED}"
fi
