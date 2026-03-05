#!/usr/bin/env bash
# =============================================================================
# scripts/05_magma.sh  —  Step V (sub-step B): MAGMA gene & gene-set enrichment
# =============================================================================
#
# Called from run.sh Step V for each stratum.
# Usage:  bash scripts/05_magma.sh <stratum>   (stratum defaults to "combined")
#
# Performs four sub-steps for each stratum:
#   5M-1. Prepare MAGMA input files (snp_loc.txt, pval.txt) via 05_magma_prep.py
#   5M-2. Annotate SNPs to genes (magma --annotate)
#   5M-3. Gene-based association test  (magma --gene-model snp-wise=mean)
#   5M-4. Gene-set analyses for Hallmark, C2, C5 (if MAGMA_GENESETS_DIR is set)
#   5M-5. Visualisation + summary tables (05_magma_visualize.py)
#
# Required config keys (from config.conf):
#   MAGMA_PATH         — path to the magma binary (auto-detected from PATH if empty)
#   MAGMA_REF_DIR      — dir containing NCBI38.gene.loc and g1000_eur.{bed,bim,fam}
#   MAGMA_GENESETS_DIR — dir with MSigDB .gmt files (Hallmark/C2/C5); leave empty
#                         to skip gene-set analyses
#   MAGMA_SAMPLE_N     — sample size N (optional; read from run_stats.tsv if absent)
#   MAGMA_GENE_WINDOW_KB — up/downstream window for SNP→gene annotation (default: 10)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

STRATUM="${1:-combined}"
STRATUM_UC=$(echo "${STRATUM}" | tr 'a-z' 'A-Z')

# ── Re-read metadata from previous steps ─────────────────────────────────────
META="${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
if [ -f "${META}" ]; then
    while IFS='=' read -r _k _v; do
        [[ "${_k}" =~ ^# ]] && continue
        [[ -z "${_k}" ]] && continue
        export "${_k}=${_v}"
    done < "${META}"
fi

# ── Resolve stratum-specific merged results ───────────────────────────────────
if [ "${STRATUM}" = "combined" ]; then
    eval "STRATUM_MERGED=\"\${MERGED_RESULTS:-}\""
else
    eval "STRATUM_MERGED=\"\${MERGED_RESULTS_${STRATUM_UC}:-}\""
fi

if [ -z "${STRATUM_MERGED:-}" ] || [ ! -f "${STRATUM_MERGED}" ]; then
    log_warn "  MAGMA [${STRATUM}]: merged results not found — skipping"
    log_warn "  (expected key MERGED_RESULTS${STRATUM:+_${STRATUM_UC}} in metadata)"
    exit 0
fi

# ── Locate MAGMA binary ───────────────────────────────────────────────────────
if [ -n "${MAGMA_PATH:-}" ] && [ -x "${MAGMA_PATH}" ]; then
    MAGMA_BIN="${MAGMA_PATH}"
elif command -v magma &>/dev/null; then
    MAGMA_BIN=$(command -v magma)
else
    log_warn "  MAGMA binary not found. Set MAGMA_PATH in config. Skipping MAGMA."
    exit 0
fi
log_info "  MAGMA binary : ${MAGMA_BIN}"

# ── Validate reference directory ─────────────────────────────────────────────
REF_DIR="${MAGMA_REF_DIR:-}"
if [ -z "${REF_DIR}" ] || [ ! -d "${REF_DIR}" ]; then
    log_warn "  MAGMA_REF_DIR is not set or does not exist. Skipping MAGMA."
    log_warn "  Set MAGMA_REF_DIR to a directory containing NCBI38.gene.loc and g1000_eur.*"
    exit 0
fi

GENE_LOC=$(find "${REF_DIR}" -name "NCBI*.gene.loc" -maxdepth 2 | sort -V | head -1)
REF_BED=$(find "${REF_DIR}"  -name "g1000_eur.bed"  -maxdepth 2 | head -1)
REF_BFILE="${REF_BED%.bed}"

if [ -z "${GENE_LOC}" ]; then
    log_warn "  NCBI*.gene.loc not found in ${REF_DIR}. Skipping MAGMA."
    exit 0
fi
if [ -z "${REF_BED}" ]; then
    log_warn "  g1000_eur.bed not found in ${REF_DIR}. Skipping MAGMA."
    exit 0
fi
log_info "  Gene-loc file: ${GENE_LOC}"
log_info "  Ref bfile    : ${REF_BFILE}"

# ── Directories ───────────────────────────────────────────────────────────────
MAGMA_DIR="${OUTPUT_DIR}/magma/${STRATUM}"
MAGMA_INPUT="${MAGMA_DIR}/input"
MAGMA_OUTPUT="${MAGMA_DIR}/output"
MAGMA_PLOTS="${MAGMA_DIR}/plots"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${MAGMA_INPUT}" "${MAGMA_OUTPUT}" "${MAGMA_PLOTS}" "${LOG_DIR}"

OUT_PREFIX="${MAGMA_OUTPUT}/${RUN_NAME}_${STRATUM}"
LOG_FILE="${LOG_DIR}/magma_${STRATUM}.log"

log_section "[MAGMA] Stratum: ${STRATUM}"

# ── 5M-1. Prepare MAGMA input files ──────────────────────────────────────────
log_info "  [5M-1] Preparing MAGMA input files"
python3 "${SCRIPT_DIR}/scripts/05_magma_prep.py" \
    --merged-results "${STRATUM_MERGED}" \
    --pval-col       "${PVAL_COL:-P_manta}" \
    --snp-loc        "${MAGMA_INPUT}/snp_loc.txt" \
    --pval-out       "${MAGMA_INPUT}/pval.txt" \
    2>&1 | tee -a "${LOG_FILE}"

N_SNPS=$(count_data_rows "${MAGMA_INPUT}/snp_loc.txt")
log_ok "  [5M-1] ${N_SNPS} SNPs with valid rsID + p-value"
if [ "${N_SNPS}" -lt 500 ]; then
    log_warn "  Only ${N_SNPS} SNPs have rsIDs — results may be sparse."
    log_warn "  Ensure Step IV rsID annotation ran successfully."
fi

# ── Resolve sample size N ─────────────────────────────────────────────────────
# Try (in order): run_stats.tsv, then MAGMA_SAMPLE_N config key
N_SAMPLES=""
STATS_FILE="${OUTPUT_DIR}/run_stats.tsv"
if [ -f "${STATS_FILE}" ]; then
    _N=$(grep "^n_samples_${STRATUM}" "${STATS_FILE}" | cut -f2 | head -1 || true)
    [ -z "${_N}" ] && _N=$(grep "^n_samples$" "${STATS_FILE}" | cut -f2 | head -1 || true)
    [ -z "${_N}" ] && _N=$(grep "^n_samples_combined" "${STATS_FILE}" | cut -f2 | head -1 || true)
    N_SAMPLES="${_N}"
fi
[ -z "${N_SAMPLES}" ] && N_SAMPLES="${MAGMA_SAMPLE_N:-}"

if [ -n "${N_SAMPLES}" ] && [[ "${N_SAMPLES}" =~ ^[0-9]+$ ]] && [ "${N_SAMPLES}" -gt 0 ]; then
    NOBS_ARG="nobs=${N_SAMPLES}"
    log_info "  Sample size  : N=${N_SAMPLES}"
else
    NOBS_ARG=""
    log_warn "  Sample size N unknown; set MAGMA_SAMPLE_N in config for accurate results"
fi

WINDOW_KB="${MAGMA_GENE_WINDOW_KB:-10}"

# ── 5M-2. Annotate SNPs → genes ──────────────────────────────────────────────
ANNOT_FILE="${OUT_PREFIX}.genes.annot"
if [ -f "${ANNOT_FILE}" ] && [ -s "${ANNOT_FILE}" ]; then
    log_info "  [5M-2] Annotation file already exists — skipping"
else
    log_info "  [5M-2] Annotating SNPs to genes (±${WINDOW_KB} kb window)"
    if [ "${DRY_RUN:-false}" = "true" ]; then
        log_info "  [DRY-RUN] Would run: magma --annotate window=${WINDOW_KB} --snp-loc ... --gene-loc ... --out ${OUT_PREFIX}"
    else
        "${MAGMA_BIN}" \
            --annotate window="${WINDOW_KB}" \
            --snp-loc  "${MAGMA_INPUT}/snp_loc.txt" \
            --gene-loc "${GENE_LOC}" \
            --out      "${OUT_PREFIX}" \
            2>&1 | tee -a "${LOG_FILE}"
        log_ok "  [5M-2] Annotation complete: ${ANNOT_FILE}"
    fi
fi

# ── 5M-3. Gene-based association test ────────────────────────────────────────
GENES_OUT="${OUT_PREFIX}.genes.out"
GENES_RAW="${OUT_PREFIX}.genes.raw"
if [ -f "${GENES_OUT}" ] && [ -s "${GENES_OUT}" ]; then
    log_info "  [5M-3] Gene-based results already exist — skipping"
else
    log_info "  [5M-3] Running gene-based test (snp-wise mean model)"
    if [ "${DRY_RUN:-false}" = "true" ]; then
        log_info "  [DRY-RUN] Would run: magma --bfile ${REF_BFILE} --gene-annot ${ANNOT_FILE} --pval ... --gene-model snp-wise=mean --out ${OUT_PREFIX}"
    else
        "${MAGMA_BIN}" \
            --bfile      "${REF_BFILE}" \
            --gene-annot "${ANNOT_FILE}" \
            --pval       "${MAGMA_INPUT}/pval.txt" use=SNP,P ${NOBS_ARG} \
            --gene-model snp-wise=mean \
            --out        "${OUT_PREFIX}" \
            2>&1 | tee -a "${LOG_FILE}"
        log_ok "  [5M-3] Gene-based test complete: ${GENES_OUT}"
    fi
fi

# ── 5M-4. Gene-set analyses ───────────────────────────────────────────────────
GENESETS_DIR="${MAGMA_GENESETS_DIR:-}"
if [ -z "${GENESETS_DIR}" ] || [ ! -d "${GENESETS_DIR}" ]; then
    log_warn "  [5M-4] MAGMA_GENESETS_DIR not set or missing — skipping gene-set analysis"
    log_warn "  Download MSigDB .gmt files (Hallmark/C2/C5) from https://www.gsea-msigdb.org"
    log_warn "  and set MAGMA_GENESETS_DIR in your config."
elif [ ! -f "${GENES_RAW}" ]; then
    log_warn "  [5M-4] genes.raw not found (gene-based test may have failed) — skipping gene-set analysis"
else
    for COLLECTION in h c2 c5; do
        # Match: h.all.*.Hs.entrez.gmt  /  c2.all.*.Hs.entrez.gmt  etc.
        GMT=$(find "${GENESETS_DIR}" -name "${COLLECTION}.all.*.Hs.entrez.gmt" | sort -V | tail -1)
        [ -z "${GMT}" ] && GMT=$(find "${GENESETS_DIR}" -name "${COLLECTION}.*.entrez.gmt" | sort -V | tail -1)
        if [ -z "${GMT}" ]; then
            log_warn "  [5M-4] No .gmt for collection '${COLLECTION}' in ${GENESETS_DIR} — skipping"
            continue
        fi
        GSA_OUT="${OUT_PREFIX}_${COLLECTION}.gsa.out"
        if [ -f "${GSA_OUT}" ] && [ -s "${GSA_OUT}" ]; then
            log_info "  [5M-4] ${COLLECTION^^} gene-set results already exist — skipping"
            continue
        fi
        log_info "  [5M-4] Gene-set analysis: ${COLLECTION^^}  ($(basename "${GMT}"))"
        if [ "${DRY_RUN:-false}" = "true" ]; then
            log_info "  [DRY-RUN] Would run: magma --gene-results ${GENES_RAW} --set-annot ${GMT} --out ${OUT_PREFIX}_${COLLECTION}"
        else
            "${MAGMA_BIN}" \
                --gene-results "${GENES_RAW}" \
                --set-annot    "${GMT}" \
                --out          "${OUT_PREFIX}_${COLLECTION}" \
                2>&1 | tee -a "${LOG_FILE}"
            log_ok "  [5M-4] ${COLLECTION^^} done: ${GSA_OUT}"
        fi
    done
fi

# ── 5M-5. Visualisation and summary tables ────────────────────────────────────
log_info "  [5M-5] Generating MAGMA plots and summary tables"
if [ "${DRY_RUN:-false}" = "false" ] && [ -f "${GENES_OUT}" ]; then
    python3 "${SCRIPT_DIR}/scripts/05_magma_visualize.py" \
        --output-prefix "${OUT_PREFIX}" \
        --plots-dir     "${MAGMA_PLOTS}" \
        --gene-loc      "${GENE_LOC}" \
        --run-name      "${RUN_NAME}_${STRATUM}" \
        2>&1 | tee -a "${LOG_DIR}/magma_visualize_${STRATUM}.log"
    log_ok "  [5M-5] Plots saved to ${MAGMA_PLOTS}/"
elif [ "${DRY_RUN:-false}" = "true" ]; then
    log_info "  [DRY-RUN] Would generate MAGMA plots"
else
    log_warn "  [5M-5] genes.out not found — skipping visualization"
fi

# ── Export paths to metadata ──────────────────────────────────────────────────
TAG="${STRATUM_UC:+_${STRATUM_UC}}"
{
    echo "MAGMA_GENES_OUT${TAG}=${GENES_OUT}"
    echo "MAGMA_PLOTS_DIR${TAG}=${MAGMA_PLOTS}"
    echo "MAGMA_OUTPUT_DIR${TAG}=${MAGMA_OUTPUT}"
} >> "${META}"

log_ok "  MAGMA complete [stratum: ${STRATUM}]"
log_info "  Gene results  : ${GENES_OUT}"
log_info "  Plots dir     : ${MAGMA_PLOTS}/"
log_info "  Log           : ${LOG_FILE}"
