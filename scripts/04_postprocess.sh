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

# Read strata from metadata
read -ra STRATA_ARR <<< "${STRATA:-combined}"
log_info "  Strata to post-process: ${STRATA_ARR[*]}"

PVAL_COL="${PVAL_COL:-P_manta}"
TOP_N="${TOP_N_SNPS:-1000}"
GENOME_BUILD="${GENOME_BUILD:-GRCh38}"
GW_THRESH="${GW_THRESH:-5e-8}"
SUG_THRESH="${SUG_THRESH:-1e-5}"
TOP_LABEL_N="${TOP_LABEL_N:-20}"
REGIONAL_WINDOW="${REGIONAL_WINDOW:-2000000}"
LOGS_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${LOGS_DIR}"

# Determine dbSNP source once — shared across all strata
if [ -n "${DBSNP_VCF:-}" ] && [ -f "${DBSNP_VCF}" ]; then
    DBSNP_SOURCE="${DBSNP_VCF}"
    log_info "  Using local dbSNP VCF: ${DBSNP_VCF}"
else
    [ "${GENOME_BUILD}" = "GRCh38" ] && \
        DBSNP_SOURCE="https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz" || \
        DBSNP_SOURCE="https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz"
    log_warn "  DBSNP_VCF not set. Falling back to remote: ${DBSNP_SOURCE}"
    log_warn "  Set DBSNP_VCF in your config for faster annotation."
fi

# =============================================================================
# Per-stratum loop: 4A top SNPs → 4B rsID annotation → 4C plots
# =============================================================================
for STRATUM in "${STRATA_ARR[@]}"; do
    log_section "[Step IV] Stratum: ${STRATUM}"

    STRATUM_UC=$(echo "${STRATUM}" | tr 'a-z' 'A-Z')
    if [ "${STRATUM}" = "combined" ]; then
        eval "STRATUM_MERGED=\"\${MERGED_RESULTS:-}\""
        STRATUM_RESULTS_DIR="${OUTPUT_DIR}/results"
        STRATUM_PLOTS_DIR="${OUTPUT_DIR}/plots"
        STRATUM_RUN_LABEL="${RUN_NAME}"
        STRATUM_STAT_SUFFIX=""
    else
        eval "STRATUM_MERGED=\"\${MERGED_RESULTS_${STRATUM_UC}:-}\""
        STRATUM_RESULTS_DIR="${OUTPUT_DIR}/results_${STRATUM}"
        STRATUM_PLOTS_DIR="${OUTPUT_DIR}/plots_${STRATUM}"
        STRATUM_RUN_LABEL="${RUN_NAME}_${STRATUM}"
        STRATUM_STAT_SUFFIX="_${STRATUM}"
    fi

    if [ -z "${STRATUM_MERGED}" ] || [ ! -f "${STRATUM_MERGED}" ]; then
        log_warn "  No merged results for stratum '${STRATUM}' — skipping (expected: ${STRATUM_MERGED:-<unset>})"
        continue
    fi
    mkdir -p "${STRATUM_RESULTS_DIR}" "${STRATUM_PLOTS_DIR}"
    log_info "  Merged results : ${STRATUM_MERGED}"

    # ── 4A. Extract top N SNPs ───────────────────────────────────────────────
    log_info "  [4A] Extracting top ${TOP_N} SNPs by ${PVAL_COL}"
    TOP_SNPS="${STRATUM_RESULTS_DIR}/top_${TOP_N}_snps.tsv"
if [ -f "${TOP_SNPS}" ] && [ -s "${TOP_SNPS}" ]; then
    log_info "  Top SNPs file already exists — skipping"
else
    if [ "${DRY_RUN}" = "false" ]; then
        # Determine column index of PVAL_COL
        HEADER=$(head -1 "${STRATUM_MERGED}")
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
        head -1 "${STRATUM_MERGED}" > "${TOP_SNPS}"
        awk -v col="${COL_IDX}" 'NR>1 && $col!="" && $col!="NA" && $col+0==$col' \
            "${STRATUM_MERGED}" | \
            sort -k"${COL_IDX},${COL_IDX}"g | \
            head -n "${TOP_N}" >> "${TOP_SNPS}"

        N_TOP=$(( $(wc -l < "${TOP_SNPS}") - 1 ))
        log_ok "  Top ${N_TOP} associations → ${TOP_SNPS}"
        record_stat "top_n_extracted${STRATUM_STAT_SUFFIX}" "${N_TOP}"
        echo "TOP_SNPS${STRATUM_UC:+_${STRATUM_UC}}=${TOP_SNPS}" >> "${META}"
    else
        log_info "  [DRY-RUN] Would extract top ${TOP_N} SNPs"
    fi
fi

# =============================================================================
# 4B. rsID annotation
# =============================================================================
log_section "[Step IV] rsID annotation (build: ${GENOME_BUILD})"

    # ── 4B. rsID annotation ──────────────────────────────────────────────────
    log_info "  [4B] rsID annotation (build: ${GENOME_BUILD})"
    ANNOTATED_TSV="${STRATUM_RESULTS_DIR}/top_${TOP_N}_snps_rsid.tsv"

if [ -f "${ANNOTATED_TSV}" ] && [ -s "${ANNOTATED_TSV}" ]; then
    log_info "  Annotated file already exists — skipping annotation"
elif [ "${DRY_RUN}" = "false" ] && [ -f "${TOP_SNPS}" ]; then
    # ── Build a minimal VCF from the top SNPs for bcftools annotate ──────────
    TMP_VCF="${STRATUM_RESULTS_DIR}/.top_snps_tmp.vcf"
    TMP_ANNOTATED="${STRATUM_RESULTS_DIR}/.top_snps_annotated.vcf"

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
    LOOKUP="${STRATUM_RESULTS_DIR}/.rsid_lookup.tsv"
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
    record_stat "rsid_annotated${STRATUM_STAT_SUFFIX}" "true"
    echo "ANNOTATED_SNPS${STRATUM_UC:+_${STRATUM_UC}}=${ANNOTATED_TSV}" >> "${META}"
    rm -f "${TMP_VCF}.gz" "${TMP_VCF}.gz.tbi" "${TMP_ANNOTATED}.gz" "${TMP_ANNOTATED}.gz.tbi" "${LOOKUP}"
else
    log_info "  [DRY-RUN] Would annotate rsIDs for stratum '${STRATUM}'"
fi

    # Use annotated file if available, else top SNPs
    PLOT_TOP="${ANNOTATED_TSV:-${TOP_SNPS}}"
    [ -f "${PLOT_TOP}" ] || PLOT_TOP="${TOP_SNPS}"

    # ── 4C. Plots ────────────────────────────────────────────────────────────
    log_info "  [4C] Generating plots for stratum '${STRATUM}'"
    if [ "${DRY_RUN}" = "false" ]; then
        Rscript "${SCRIPT_DIR}/scripts/04_plot_manhattan.R" \
            --input        "${STRATUM_MERGED}" \
            --top-file     "${PLOT_TOP}" \
            --output-dir   "${STRATUM_PLOTS_DIR}" \
            --run-name     "${STRATUM_RUN_LABEL}" \
            --pval-col     "${PVAL_COL}" \
            --genome-build "${GENOME_BUILD}" \
            --gw-thresh    "${GW_THRESH}" \
            --sug-thresh   "${SUG_THRESH}" \
            --top-label-n  "${TOP_LABEL_N}" \
            --regional-window "${REGIONAL_WINDOW}" \
            2>&1 | tee -a "${LOGS_DIR}/plotting_${STRATUM}.log"

        for PLOT_TYPE in manhattan qq regional; do
            PNG="${STRATUM_PLOTS_DIR}/${PLOT_TYPE}_${STRATUM_RUN_LABEL}.png"
            [ -f "${PNG}" ] \
                && log_ok "  ${PLOT_TYPE^}: ${PNG}" \
                || log_warn "  ${PLOT_TYPE^} not found — check ${LOGS_DIR}/plotting_${STRATUM}.log"
        done

        echo "PLOT_MANHATTAN${STRATUM_UC:+_${STRATUM_UC}}=${STRATUM_PLOTS_DIR}/manhattan_${STRATUM_RUN_LABEL}.png" >> "${META}"
        echo "PLOT_QQ${STRATUM_UC:+_${STRATUM_UC}}=${STRATUM_PLOTS_DIR}/qq_${STRATUM_RUN_LABEL}.png"             >> "${META}"
        echo "PLOT_REGIONAL${STRATUM_UC:+_${STRATUM_UC}}=${STRATUM_PLOTS_DIR}/regional_${STRATUM_RUN_LABEL}.png" >> "${META}"
        echo "PLOTS_DIR${STRATUM_UC:+_${STRATUM_UC}}=${STRATUM_PLOTS_DIR}"                                       >> "${META}"
    else
        log_info "  [DRY-RUN] Would call 04_plot_manhattan.R for stratum '${STRATUM}'"
    fi

done  # end strata loop

log_milestone "Step IV complete — top SNPs, rsID annotation and plots done for all strata"
