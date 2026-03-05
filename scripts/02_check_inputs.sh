#!/usr/bin/env bash
# =============================================================================
# scripts/02_check_inputs.sh  —  Step II: Input validation and ID harmonisation
# =============================================================================
# Validates the three input files (VCF, phenotype, covariate) and harmonises
# the sample IDs so that only the intersection of all three is retained.
#
# Outputs written to ${OUTPUT_DIR}/inputs/:
#   samples_vcf.txt            — sample IDs from VCF
#   samples_pheno.txt          — sample IDs from phenotype file
#   samples_cov.txt            — sample IDs from covariate file
#   samples_intersection.txt   — overlapping IDs (used downstream)
#   phenotype_filtered.tsv     — phenotype file filtered to intersection
#   covariate_filtered.tsv     — covariate file filtered to intersection
#   input_validation_report.txt — detailed text report
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

INPUTS_DIR="${OUTPUT_DIR}/inputs"
mkdir -p "${INPUTS_DIR}"

REPORT="${INPUTS_DIR}/input_validation_report.txt"
: > "${REPORT}"

_report() { echo "$@" | tee -a "${REPORT}"; }

_report "============================================================"
_report "  Input Validation Report — ${RUN_NAME}"
_report "  Generated: $(date '+%Y-%m-%d %H:%M:%S')"
_report "============================================================"

# =============================================================================
# 2A. Check input files exist
# =============================================================================
log_section "[Step II] File existence checks"

require_file "${VCF_FILE}"        "VCF file"
require_file "${PHENOTYPE_FILE}"  "Phenotype file"
require_file "${COVARIATE_FILE}"  "Covariate file"

_report ""
_report "Input files"
_report "  VCF         : ${VCF_FILE}"
_report "  Phenotype   : ${PHENOTYPE_FILE}"
_report "  Covariate   : ${COVARIATE_FILE}"

# =============================================================================
# 2B. VCF file validation
# =============================================================================
log_section "[Step II] VCF file validation"

# Check index
require_index "${VCF_FILE}"

# Extract sample IDs from VCF
log_info "  Extracting sample IDs from VCF..."
bcftools query -l "${VCF_FILE}" > "${INPUTS_DIR}/samples_vcf.txt"
N_VCF=$(wc -l < "${INPUTS_DIR}/samples_vcf.txt")
log_ok "  VCF samples: $(fmt_num ${N_VCF})"

# Check VCF has at least a few SNPs by counting chromosomes present
log_info "  Checking chromosomes in VCF..."
CHR_LIST=$(bcftools view -H "${VCF_FILE}" 2>/dev/null | awk '{print $1}' | sort -u | head -50 || \
           bcftools index -s "${VCF_FILE}" | awk '{print $1}' | head -50)
log_info "  Chromosomes present (first 10): $(echo "${CHR_LIST}" | head -10 | tr '\n' ' ')"

# Detect chromosome naming: chr-prefix vs plain numbers
if echo "${CHR_LIST}" | grep -q '^chr'; then
    VCF_CHR_PREFIX="chr"
    log_info "  Chromosome naming: chr-prefixed (e.g. chr1, chr2...)"
else
    VCF_CHR_PREFIX=""
    log_info "  Chromosome naming: plain numbers (e.g. 1, 2...)"
fi
echo "VCF_CHR_PREFIX=${VCF_CHR_PREFIX}" >> "${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"

_report ""
_report "VCF file"
_report "  Path         : ${VCF_FILE}"
_report "  Samples      : ${N_VCF}"
_report "  Chr prefix   : '${VCF_CHR_PREFIX}'"
_report "  Chromosomes  : $(echo "${CHR_LIST}" | wc -w | tr -d ' ')"

record_stat "n_samples_vcf" "${N_VCF}"

# =============================================================================
# 2C. Phenotype file validation
# =============================================================================
log_section "[Step II] Phenotype file validation"

PHENO_HEADER=$(get_header "${PHENOTYPE_FILE}")
N_PHENO_ROWS=$(count_data_rows "${PHENOTYPE_FILE}")
N_PHENO_COLS=$(echo "${PHENO_HEADER}" | awk -F'\t' '{print NF}')
N_PHENOTYPES=$(( N_PHENO_COLS - 1 ))

log_info "  Phenotype rows (data): $(fmt_num ${N_PHENO_ROWS})"
log_info "  Phenotype columns: ${N_PHENO_COLS}  (ID + ${N_PHENOTYPES} phenotypes)"
log_info "  First few column names: $(echo "${PHENO_HEADER}" | cut -f1-6 | tr '\t' '|')"

# Detect ID column: must be first column
PHENO_ID_COL=$(echo "${PHENO_HEADER}" | awk -F'\t' '{print $1}')
log_info "  ID column: '${PHENO_ID_COL}'"

# Extract sample IDs
if [[ "${PHENOTYPE_FILE}" == *.gz ]]; then
    zcat "${PHENOTYPE_FILE}" | tail -n +2 | awk -F'\t' '{print $1}' > "${INPUTS_DIR}/samples_pheno.txt"
else
    tail -n +2 "${PHENOTYPE_FILE}" | awk -F'\t' '{print $1}' > "${INPUTS_DIR}/samples_pheno.txt"
fi
N_PHENO_IDS=$(wc -l < "${INPUTS_DIR}/samples_pheno.txt")

# Check for duplicate IDs in phenotype file
N_PHENO_UNIQ=$(sort -u "${INPUTS_DIR}/samples_pheno.txt" | wc -l)
if [ "${N_PHENO_IDS}" -ne "${N_PHENO_UNIQ}" ]; then
    log_warn "  Duplicate sample IDs in phenotype file! (${N_PHENO_IDS} rows, ${N_PHENO_UNIQ} unique)"
    _report "  WARNING: Duplicate IDs in phenotype file"
else
    log_ok "  No duplicate IDs in phenotype file."
fi

# Check for missing values
if [[ "${PHENOTYPE_FILE}" == *.gz ]]; then
    N_NA=$(zcat "${PHENOTYPE_FILE}" | tail -n +2 | grep -cE '\t(NA|na|nan|NaN|\.)(\t|$)' || echo 0)
else
    N_NA=$(tail -n +2 "${PHENOTYPE_FILE}" | grep -cE '\t(NA|na|nan|NaN|\.)(\t|$)' || echo 0)
fi
if [ "${N_NA}" -gt 0 ]; then
    log_warn "  ${N_NA} rows contain NA/missing values in phenotype file."
    _report "  WARNING: ${N_NA} rows with NA in phenotype file"
else
    log_ok "  No NA values detected in phenotype file."
fi

_report ""
_report "Phenotype file"
_report "  Path          : ${PHENOTYPE_FILE}"
_report "  Data rows     : ${N_PHENO_ROWS}"
_report "  Phenotypes    : ${N_PHENOTYPES}"
_report "  ID column     : ${PHENO_ID_COL}"
_report "  Unique IDs    : ${N_PHENO_UNIQ}"

record_stat "n_phenotypes"    "${N_PHENOTYPES}"
record_stat "n_samples_pheno" "${N_PHENO_IDS}"

# =============================================================================
# 2D. Covariate file validation
# =============================================================================
log_section "[Step II] Covariate file validation"

COV_HEADER=$(get_header "${COVARIATE_FILE}")
N_COV_ROWS=$(count_data_rows "${COVARIATE_FILE}")
N_COV_COLS=$(echo "${COV_HEADER}" | awk -F'\t' '{print NF}')
N_COVARIATES=$(( N_COV_COLS - 1 ))

log_info "  Covariate rows (data): $(fmt_num ${N_COV_ROWS})"
log_info "  Covariate columns: ${N_COV_COLS}  (ID + ${N_COVARIATES} covariates)"
log_info "  Covariate names: $(echo "${COV_HEADER}" | cut -f2- | tr '\t' ', ')"

# Extract sample IDs
if [[ "${COVARIATE_FILE}" == *.gz ]]; then
    zcat "${COVARIATE_FILE}" | tail -n +2 | awk -F'\t' '{print $1}' > "${INPUTS_DIR}/samples_cov.txt"
else
    tail -n +2 "${COVARIATE_FILE}" | awk -F'\t' '{print $1}' > "${INPUTS_DIR}/samples_cov.txt"
fi
N_COV_IDS=$(wc -l < "${INPUTS_DIR}/samples_cov.txt")

N_COV_UNIQ=$(sort -u "${INPUTS_DIR}/samples_cov.txt" | wc -l)
if [ "${N_COV_IDS}" -ne "${N_COV_UNIQ}" ]; then
    log_warn "  Duplicate sample IDs in covariate file! (${N_COV_IDS} rows, ${N_COV_UNIQ} unique)"
else
    log_ok "  No duplicate IDs in covariate file."
fi

_report ""
_report "Covariate file"
_report "  Path          : ${COVARIATE_FILE}"
_report "  Data rows     : ${N_COV_ROWS}"
_report "  Covariates    : ${N_COVARIATES}"
_report "  Unique IDs    : ${N_COV_UNIQ}"

record_stat "n_covariates"    "${N_COVARIATES}"
record_stat "n_samples_cov"   "${N_COV_IDS}"

# =============================================================================
# 2E. Compute intersection of sample IDs
# =============================================================================
log_section "[Step II] Sample ID intersection (VCF ∩ phenotype ∩ covariate)"

# Sort all three lists
sort "${INPUTS_DIR}/samples_vcf.txt"   > "${INPUTS_DIR}/samples_vcf_sorted.txt"
sort "${INPUTS_DIR}/samples_pheno.txt" > "${INPUTS_DIR}/samples_pheno_sorted.txt"
sort "${INPUTS_DIR}/samples_cov.txt"   > "${INPUTS_DIR}/samples_cov_sorted.txt"

# Two-way comm operations to get full intersection
comm -12 "${INPUTS_DIR}/samples_vcf_sorted.txt" \
         "${INPUTS_DIR}/samples_pheno_sorted.txt" \
         > "${INPUTS_DIR}/samples_vcf_pheno_common.txt"

comm -12 "${INPUTS_DIR}/samples_vcf_pheno_common.txt" \
         "${INPUTS_DIR}/samples_cov_sorted.txt" \
         > "${INPUTS_DIR}/samples_intersection.txt"

N_INTERSECT=$(wc -l < "${INPUTS_DIR}/samples_intersection.txt")

log_info "  VCF samples           : $(fmt_num ${N_VCF})"
log_info "  Phenotype samples     : $(fmt_num ${N_PHENO_IDS})"
log_info "  Covariate samples     : $(fmt_num ${N_COV_IDS})"
log_ok   "  Intersection (all 3)  : $(fmt_num ${N_INTERSECT})"

# Check what will be excluded
N_ONLY_VCF=$(comm -23 "${INPUTS_DIR}/samples_vcf_sorted.txt" \
                       "${INPUTS_DIR}/samples_intersection.txt" | wc -l)
N_ONLY_PHENO=$(comm -23 "${INPUTS_DIR}/samples_pheno_sorted.txt" \
                         "${INPUTS_DIR}/samples_intersection.txt" | wc -l)
N_ONLY_COV=$(comm -23 "${INPUTS_DIR}/samples_cov_sorted.txt" \
                       "${INPUTS_DIR}/samples_intersection.txt" | wc -l)

if [ "${N_ONLY_VCF}" -gt 0 ]; then
    log_warn "  ${N_ONLY_VCF} VCF samples not in phenotype/covariate files (will be excluded at bcftools step)"
fi
if [ "${N_ONLY_PHENO}" -gt 0 ]; then
    log_warn "  ${N_ONLY_PHENO} phenotype samples not in VCF/covariate (will be dropped from filtered files)"
fi
if [ "${N_ONLY_COV}" -gt 0 ]; then
    log_warn "  ${N_ONLY_COV} covariate samples not in VCF/phenotype (will be dropped from filtered files)"
fi

if [ "${N_INTERSECT}" -lt 50 ]; then
    die "Intersection contains only ${N_INTERSECT} samples — too few for GWAS. Check sample ID format across files."
fi

_report ""
_report "Sample ID overlap"
_report "  VCF samples       : ${N_VCF}"
_report "  Phenotype samples : ${N_PHENO_IDS}"
_report "  Covariate samples : ${N_COV_IDS}"
_report "  INTERSECTION      : ${N_INTERSECT}"
_report "  Excluded (VCF only)    : ${N_ONLY_VCF}"
_report "  Excluded (pheno only)  : ${N_ONLY_PHENO}"
_report "  Excluded (cov only)    : ${N_ONLY_COV}"

record_stat "n_samples_intersection" "${N_INTERSECT}"

# =============================================================================
# 2F. Filter phenotype and covariate files to intersection
# =============================================================================
log_section "[Step II] Filtering input files to intersection"

python3 "${SCRIPT_DIR}/scripts/02_filter_individuals.py" \
    --pheno-in       "${PHENOTYPE_FILE}" \
    --cov-in         "${COVARIATE_FILE}" \
    --samples        "${INPUTS_DIR}/samples_intersection.txt" \
    --pheno-out      "${INPUTS_DIR}/phenotype_filtered.tsv" \
    --cov-out        "${INPUTS_DIR}/covariate_filtered.tsv" \
    --report         "${INPUTS_DIR}/filter_report.txt" \
    --verbose

require_file "${INPUTS_DIR}/phenotype_filtered.tsv" "Filtered phenotype file"
require_file "${INPUTS_DIR}/covariate_filtered.tsv" "Filtered covariate file"

N_PHENO_FILTERED=$(count_data_rows "${INPUTS_DIR}/phenotype_filtered.tsv")
N_COV_FILTERED=$(count_data_rows "${INPUTS_DIR}/covariate_filtered.tsv")
log_ok "  Filtered phenotype: ${N_PHENO_FILTERED} samples"
log_ok "  Filtered covariate: ${N_COV_FILTERED} samples"

if [ "${N_PHENO_FILTERED}" -ne "${N_INTERSECT}" ] || [ "${N_COV_FILTERED}" -ne "${N_INTERSECT}" ]; then
    die "Post-filter sample counts do not match intersection (expected ${N_INTERSECT}). Something went wrong."
fi

# Export paths for use by downstream steps
echo "PHENOTYPE_FILTERED=${INPUTS_DIR}/phenotype_filtered.tsv" >> "${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
echo "COVARIATE_FILTERED=${INPUTS_DIR}/covariate_filtered.tsv" >> "${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
echo "SAMPLES_INTERSECTION=${INPUTS_DIR}/samples_intersection.txt" >> "${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"
echo "VCF_CHR_PREFIX=${VCF_CHR_PREFIX}" >> "${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"

_report ""
_report "Filtered files (used in pipeline)"
_report "  Phenotype   : ${INPUTS_DIR}/phenotype_filtered.tsv  (${N_PHENO_FILTERED} samples)"
_report "  Covariate   : ${INPUTS_DIR}/covariate_filtered.tsv  (${N_COV_FILTERED} samples)"
_report "  Sample list : ${INPUTS_DIR}/samples_intersection.txt"

_report ""
_report "============================================================"
_report "  VALIDATION COMPLETE — all checks passed"
_report "  Analysis-ready sample size: ${N_INTERSECT}"
_report "============================================================"

log_ok "Step II complete — ${N_INTERSECT} samples harmonised across all three input files."
