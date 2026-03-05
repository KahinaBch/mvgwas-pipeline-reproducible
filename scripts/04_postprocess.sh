#!/usr/bin/env bash
# =============================================================================
# scripts/04_postprocess.sh  —  Step IV: top SNPs, rsID annotation, plots
# =============================================================================
set -euo pipefail

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

require_file "${MERGED_RESULTS:-}" "Merged GWAS results"

RESULTS_DIR="${OUTPUT_DIR}/results"
QC_DIR="${OUTPUT_DIR}/qc"
PLOTS_DIR="${OUTPUT_DIR}/plots"
LOGS_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${RESULTS_DIR}" "${QC_DIR}" "${PLOTS_DIR}"

PVAL_COL="${PVAL_COL:-P_manta}"
TOP_N="${TOP_N_SNPS:-1000}"
GENOME_BUILD="${GENOME_BUILD:-GRCh38}"
GW_THRESH="${GW_THRESH:-5e-8}"
SUG_THRESH="${SUG_THRESH:-1e-5}"
TOP_LABEL_N="${TOP_LABEL_N:-20}"
REGIONAL_WINDOW="${REGIONAL_WINDOW:-2000000}"

# =============================================================================
# 4A. Extract top N SNPs
# =============================================================================
log_section "[Step IV] Extracting top ${TOP_N} SNPs by ${PVAL_COL}"

TOP_SNPS="${RESULTS_DIR}/top_${TOP_N}_snps.tsv"
if [ -f "${TOP_SNPS}" ] && [ -s "${TOP_SNPS}" ]; then
    log_info "  Top SNPs file already exists — skipping"
else
    if [ "${DRY_RUN}" = "false" ]; then
        # Determine column index of PVAL_COL
        HEADER=$(head -1 "${MERGED_RESULTS}")
        COL_IDX=$(echo "${HEADER}" | tr '\t' '\n' | grep -n "^${PVAL_COL}$" | cut -d: -f1)
        if [ -z "${COL_IDX}" ]; then
            log_warn "  Column '${PVAL_COL}' not found — trying case-insensitive match"
            COL_IDX=$(echo "${HEADER}" | tr '\t' '\n' | grep -in "${PVAL_COL}" | head -1 | cut -d: -f1)
        fi
        if [ -z "${COL_IDX}" ]; then
            die "Cannot locate p-value column '${PVAL_COL}' in header: ${HEADER}"
        fi
        log_info "  P-value column '${PVAL_COL}' is column ${COL_IDX}"

        # Sort numerically ascending by p-value, take header + top N rows
        head -1 "${MERGED_RESULTS}" > "${TOP_SNPS}"
        awk -v col="${COL_IDX}" 'NR>1 && $col!="" && $col!="NA" && $col+0==$col' \
            "${MERGED_RESULTS}" | \
            sort -k"${COL_IDX},${COL_IDX}"g | \
            head -n "${TOP_N}" >> "${TOP_SNPS}"

        N_TOP=$(( $(wc -l < "${TOP_SNPS}") - 1 ))
        log_ok "  Top ${N_TOP} associations → ${TOP_SNPS}"
        record_stat "top_n_extracted" "${N_TOP}"
        echo "TOP_SNPS=${TOP_SNPS}" >> "${META}"
    else
        log_info "  [DRY-RUN] Would extract top ${TOP_N} SNPs"
    fi
fi

# =============================================================================
# 4B. rsID annotation
# =============================================================================
log_section "[Step IV] rsID annotation (build: ${GENOME_BUILD})"

ANNOTATED_TSV="${RESULTS_DIR}/top_${TOP_N}_snps_rsid.tsv"

# Determine dbSNP VCF source
if [ -n "${DBSNP_VCF:-}" ] && [ -f "${DBSNP_VCF}" ]; then
    log_info "  Using local dbSNP VCF: ${DBSNP_VCF}"
    DBSNP_SOURCE="${DBSNP_VCF}"
else
    # Remote URL based on genome build
    if [ "${GENOME_BUILD}" = "GRCh38" ]; then
        DBSNP_URL="https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
    else
        DBSNP_URL="https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"
    fi
    log_warn "  DBSNP_VCF not set. Falling back to remote: ${DBSNP_URL}"
    log_warn "  This may be very slow. Set DBSNP_VCF in your config for faster annotation."
    DBSNP_SOURCE="${DBSNP_URL}"
fi

if [ -f "${ANNOTATED_TSV}" ] && [ -s "${ANNOTATED_TSV}" ]; then
    log_info "  Annotated file already exists — skipping annotation"
elif [ "${DRY_RUN}" = "false" ]; then
    # ── Build a minimal VCF from the top SNPs for bcftools annotate ──────────
    TMP_VCF="${RESULTS_DIR}/.top_snps_tmp.vcf"
    TMP_ANNOTATED="${RESULTS_DIR}/.top_snps_annotated.vcf"

    HEADER=$(head -1 "${TOP_SNPS}")
    # Try MANTA-style columns (CHR/POS) or PLINK-style (#CHROM/POS/ID/REF/ALT)
    CHR_COL=$(echo "${HEADER}" | tr '\t' '\n' | grep -n "^#\?CHROM\|^CHR$" | head -1 | cut -d: -f1)
    POS_COL=$(echo "${HEADER}" | tr '\t' '\n' | grep -in "^POS$" | head -1 | cut -d: -f1)
    REF_COL=$(echo "${HEADER}" | tr '\t' '\n' | grep -in "^REF$" | head -1 | cut -d: -f1)
    ALT_COL=$(echo "${HEADER}" | tr '\t' '\n' | grep -in "^ALT$" | head -1 | cut -d: -f1)

    log_info "  CHR col: ${CHR_COL}, POS col: ${POS_COL}, REF col: ${REF_COL}, ALT col: ${ALT_COL}"

    # Create minimal VCF (CHROM POS ID REF ALT QUAL FILTER INFO)
    {
        printf "##fileformat=VCFv4.1\n"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        awk -v chr="${CHR_COL}" -v pos="${POS_COL}" -v ref="${REF_COL:-0}" \
            -v alt="${ALT_COL:-0}" 'NR>1 {
            c = $chr; sub(/^chr/, "", c)
            r = (ref>0) ? $ref : "."
            a = (alt>0) ? $alt : "."
            printf "%s\t%s\t.\t%s\t%s\t.\t.\t.\n", c, $pos, r, a
        }' "${TOP_SNPS}"
    } > "${TMP_VCF}"

    bgzip -f "${TMP_VCF}"
    tabix -p vcf "${TMP_VCF}.gz"

    log_info "  Running bcftools annotate..."
    bcftools annotate \
        -a "${DBSNP_SOURCE}" \
        -c ID \
        --regions-file <(awk 'NR>1{print $1"\t"$2"\t"$2}' "${TOP_SNPS}" | sort -k1,1 -k2,2n) \
        "${TMP_VCF}.gz" \
        -Oz -o "${TMP_ANNOTATED}.gz" 2>> "${LOGS_DIR}/rsid_annotation.log"
    tabix -p vcf "${TMP_ANNOTATED}.gz"

    # ── Build a CHR:POS → rsID lookup table ──────────────────────────────────
    LOOKUP="${RESULTS_DIR}/.rsid_lookup.tsv"
    bcftools query -f '%CHROM\t%POS\t%ID\n' "${TMP_ANNOTATED}.gz" > "${LOOKUP}"

    # ── Join rsIDs back to top SNPs table ────────────────────────────────────
    log_info "  Joining rsIDs back to association table..."
    python3 - <<PYEOF
import csv, sys

lookup = {}
with open("${LOOKUP}") as fh:
    for line in fh:
        chrom, pos, rsid = line.strip().split('\t')
        chrom = chrom.lstrip('chr')
        lookup[(chrom, pos)] = rsid if rsid != '.' else ''

with open("${TOP_SNPS}") as fin, open("${ANNOTATED_TSV}", 'w') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    # Determine CHR/POS column names
    chr_col = next((c for c in reader.fieldnames if c.upper() in ('CHR','#CHROM','CHROM')), None)
    pos_col = next((c for c in reader.fieldnames if c.upper() == 'POS'), None)
    if not chr_col or not pos_col:
        print(f"ERROR: could not find CHR/POS columns in {reader.fieldnames}", file=sys.stderr)
        sys.exit(1)

    fieldnames = ['rsID'] + reader.fieldnames
    writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for row in reader:
        chrom = row[chr_col].lstrip('chr')
        rsid = lookup.get((chrom, row[pos_col]), '')
        writer.writerow({'rsID': rsid, **row})

print(f"Annotated {sum(1 for _ in open('${ANNOTATED_TSV}'))-1} SNPs.")
PYEOF

    log_ok "  rsID annotation complete → ${ANNOTATED_TSV}"
    record_stat "rsid_annotated" "true"
    echo "ANNOTATED_SNPS=${ANNOTATED_TSV}" >> "${META}"

    # Cleanup temp files
    rm -f "${TMP_VCF}.gz" "${TMP_VCF}.gz.tbi" "${TMP_ANNOTATED}.gz" "${TMP_ANNOTATED}.gz.tbi" "${LOOKUP}"
else
    log_info "  [DRY-RUN] Would annotate rsIDs using ${DBSNP_SOURCE}"
fi

# Use annotated file if available, else top SNPs
PLOT_TOP="${ANNOTATED_TSV:-${TOP_SNPS}}"
[ -f "${PLOT_TOP}" ] || PLOT_TOP="${TOP_SNPS}"

# =============================================================================
# 4C. Generate plots: Manhattan, QQ, Regional Manhattan
# =============================================================================
log_section "[Step IV] Generating plots"

if [ "${DRY_RUN}" = "false" ]; then
    log_info "  Running Manhattan / QQ / Regional plot script..."
    Rscript "${SCRIPT_DIR}/scripts/04_plot_manhattan.R" \
        --input "${MERGED_RESULTS}" \
        --top-file "${PLOT_TOP}" \
        --output-dir "${PLOTS_DIR}" \
        --run-name "${RUN_NAME}" \
        --pval-col "${PVAL_COL}" \
        --genome-build "${GENOME_BUILD}" \
        --gw-thresh "${GW_THRESH}" \
        --sug-thresh "${SUG_THRESH}" \
        --top-label-n "${TOP_LABEL_N}" \
        --regional-window "${REGIONAL_WINDOW}" \
        2>&1 | tee -a "${LOGS_DIR}/plotting.log"

    # Verify outputs
    for PLOT_TYPE in manhattan qq regional; do
        EXPECTED_PNG="${PLOTS_DIR}/${PLOT_TYPE}_${RUN_NAME}.png"
        if [ -f "${EXPECTED_PNG}" ]; then
            log_ok "  ${PLOT_TYPE^} plot: ${EXPECTED_PNG}"
        else
            log_warn "  ${PLOT_TYPE^} plot not found — check ${LOGS_DIR}/plotting.log"
        fi
    done

    # Export plot paths to metadata
    echo "PLOT_MANHATTAN=${PLOTS_DIR}/manhattan_${RUN_NAME}.png" >> "${META}"
    echo "PLOT_QQ=${PLOTS_DIR}/qq_${RUN_NAME}.png" >> "${META}"
    echo "PLOT_REGIONAL=${PLOTS_DIR}/regional_${RUN_NAME}.png" >> "${META}"
    echo "PLOTS_DIR=${PLOTS_DIR}" >> "${META}"
else
    log_info "  [DRY-RUN] Would call 04_plot_manhattan.R"
fi

log_milestone "Step IV complete — top SNPs, rsID annotation and plots done"
