#!/usr/bin/env bash
# =============================================================================
# scripts/01_setup_test.sh  —  Step I: Environment setup and pipeline self-test
# =============================================================================
# Verifies all required tools are installed and accessible, checks Java/Nextflow
# version compatibility, then optionally runs the bundled mvgwas-nf test dataset.
#
# Called by run.sh. Expects the following env vars to be set (via config + run.sh):
#   PIPELINE_DIR, CONTAINER, NXF_VER, OUTPUT_DIR, VERBOSE, DRY_RUN
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/lib/logger.sh"
source "${SCRIPT_DIR}/lib/utils.sh"

RUN_TEST=true
if [[ "${1:-}" == "--no-test" ]]; then
    RUN_TEST=false
fi

LOG_SECTION_PREFIX="[Step I]"

# =============================================================================
# 1A. Tool availability checks
# =============================================================================
log_section "${LOG_SECTION_PREFIX} Tool availability checks"

REQUIRED_TOOLS=(java nextflow bcftools bgzip tabix Rscript python3)
OPTIONAL_TOOLS=(singularity docker)

ALL_OK=true
for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "${tool}" &>/dev/null; then
        VER=$(${tool} --version 2>&1 | head -1 || echo "(version unavailable)")
        log_ok "  ${tool}: $(command -v "${tool}")  [${VER}]"
    else
        log_error "  MISSING required tool: ${tool}"
        ALL_OK=false
    fi
done

for tool in "${OPTIONAL_TOOLS[@]}"; do
    if command -v "${tool}" &>/dev/null; then
        log_ok "  ${tool}: $(command -v "${tool}")  (optional)"
    else
        log_warn "  ${tool} not found (optional; needed if CONTAINER=${tool})"
    fi
done

if [ "${ALL_OK}" = "false" ]; then
    die "One or more required tools are missing. Install them and re-run."
fi

record_stat "step1_tools_ok" "true"

# =============================================================================
# 1B. Java version check
# =============================================================================
log_section "${LOG_SECTION_PREFIX} Java version check"

JAVA_VER=$(java -version 2>&1 | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "0.0.0")
JAVA_MAJOR=$(echo "${JAVA_VER}" | cut -d. -f1)
# Java 8 reports as 1.8, newer versions as 11, 17, 21…
if [ "${JAVA_MAJOR}" = "1" ]; then
    JAVA_MAJOR=$(echo "${JAVA_VER}" | cut -d. -f2)
fi

log_info "  Java major version: ${JAVA_MAJOR}  (full: ${JAVA_VER})"

if [ "${JAVA_MAJOR}" -ge 8 ]; then
    log_ok "  Java version is compatible (>= 8)."
else
    die "Java 8+ is required. Detected: ${JAVA_VER}"
fi

if [ "${JAVA_MAJOR}" -eq 8 ]; then
    log_warn "  Java 8 detected. Nextflow compatibility may be limited — forcing NXF_VER=22.04.0"
    export NXF_VER="${NXF_VER:-22.04.0}"
fi

record_stat "java_version" "${JAVA_VER}"

# =============================================================================
# 1C. Nextflow version check
# =============================================================================
log_section "${LOG_SECTION_PREFIX} Nextflow version check"

NXF_CMD="nextflow"
if [ -n "${NXF_VER:-}" ]; then
    log_info "  Using Nextflow version: ${NXF_VER} (set via NXF_VER)"
    export NXF_VER
fi

NXF_ACTUAL_VER=$(${NXF_CMD} -version 2>&1 | grep -oE 'version [0-9]+\.[0-9]+\.[0-9]+' | head -1 | awk '{print $2}' || echo "unknown")
log_ok "  Nextflow version: ${NXF_ACTUAL_VER}"
record_stat "nextflow_version" "${NXF_ACTUAL_VER}"

# =============================================================================
# 1D. Pipeline directory check
# =============================================================================
log_section "${LOG_SECTION_PREFIX} Pipeline directory check"

if [ ! -f "${PIPELINE_DIR}/mvgwas.nf" ]; then
    die "mvgwas.nf not found in PIPELINE_DIR=${PIPELINE_DIR}. Check your config."
fi
log_ok "  Found: ${PIPELINE_DIR}/mvgwas.nf"

# Check which branch is active
if git -C "${PIPELINE_DIR}" rev-parse --abbrev-ref HEAD &>/dev/null; then
    NF_BRANCH=$(git -C "${PIPELINE_DIR}" rev-parse --abbrev-ref HEAD)
    NF_COMMIT=$(git -C "${PIPELINE_DIR}" rev-parse --short HEAD)
    log_info "  Branch : ${NF_BRANCH}  (commit: ${NF_COMMIT})"
    if [ "${NF_BRANCH}" != "${PIPELINE_BRANCH:-kb/dsl2-conversion}" ]; then
        log_warn "  Expected branch '${PIPELINE_BRANCH:-kb/dsl2-conversion}' but found '${NF_BRANCH}'. Verify this is correct."
    else
        log_ok "  Correct branch: ${NF_BRANCH}"
    fi
    record_stat "pipeline_branch" "${NF_BRANCH}"
    record_stat "pipeline_commit" "${NF_COMMIT}"
fi

# =============================================================================
# 1E. R package checks
# =============================================================================
log_section "${LOG_SECTION_PREFIX} R package checks"

R_PACKAGES=(data.table ggplot2 ggrepel scales cowplot jsonlite optparse)
MISSING_R=()
for pkg in "${R_PACKAGES[@]}"; do
    if Rscript -e "if (!requireNamespace('${pkg}', quietly=TRUE)) quit(status=1)" 2>/dev/null; then
        log_ok "  R package ${pkg}: installed"
    else
        log_warn "  R package ${pkg}: NOT FOUND — will attempt install"
        MISSING_R+=("${pkg}")
    fi
done

if [ ${#MISSING_R[@]} -gt 0 ]; then
    log_info "  Installing missing R packages: ${MISSING_R[*]}"
    if [ "${DRY_RUN}" = "false" ]; then
        Rscript -e "install.packages(c($(printf "'%s'," "${MISSING_R[@]}" | sed 's/,$//')))" \
            --no-save --quiet \
        || log_warn "  Some R packages could not be auto-installed. Install manually with install.packages()."
    fi
fi

# =============================================================================
# 1F. Python package checks
# =============================================================================
log_section "${LOG_SECTION_PREFIX} Python package checks"

PYTHON_PACKAGES=(pandas numpy scipy matplotlib jinja2)
MISSING_PY=()
for pkg in "${PYTHON_PACKAGES[@]}"; do
    if python3 -c "import ${pkg}" 2>/dev/null; then
        VER=$(python3 -c "import ${pkg}; print(getattr(${pkg},'__version__','?'))" 2>/dev/null || echo "?")
        log_ok "  Python ${pkg}: ${VER}"
    else
        log_warn "  Python ${pkg}: NOT FOUND"
        MISSING_PY+=("${pkg}")
    fi
done

if [ ${#MISSING_PY[@]} -gt 0 ]; then
    log_info "  Installing missing Python packages: ${MISSING_PY[*]}"
    if [ "${DRY_RUN}" = "false" ]; then
        python3 -m pip install --quiet "${MISSING_PY[@]}" \
        || log_warn "  Some Python packages could not be auto-installed."
    fi
fi

# =============================================================================
# 1G. Container runtime check
# =============================================================================
log_section "${LOG_SECTION_PREFIX} Container runtime check"

case "${CONTAINER:-singularity}" in
    singularity)
        if command -v singularity &>/dev/null; then
            SIF_VER=$(singularity version 2>&1 | head -1)
            log_ok "  Singularity: ${SIF_VER}"
            NXF_CONTAINER_FLAG="-with-singularity"
        else
            log_warn "  Singularity not found. Will try Docker."
            NXF_CONTAINER_FLAG="-with-docker"
        fi
        ;;
    docker)
        if command -v docker &>/dev/null; then
            DCK_VER=$(docker --version 2>&1 | head -1)
            log_ok "  Docker: ${DCK_VER}"
            NXF_CONTAINER_FLAG="-with-docker"
        else
            die "CONTAINER=docker but docker is not installed."
        fi
        ;;
    none)
        log_warn "  CONTAINER=none — running without container isolation."
        NXF_CONTAINER_FLAG=""
        ;;
    *)
        die "Unknown CONTAINER value: ${CONTAINER}. Use singularity, docker, or none."
        ;;
esac

export NXF_CONTAINER_FLAG
echo "NXF_CONTAINER_FLAG=${NXF_CONTAINER_FLAG}" >> "${OUTPUT_DIR}/${RUN_NAME}_run_metadata.txt"

# =============================================================================
# 1H. Pipeline self-test (optional)
# =============================================================================
if [ "${RUN_TEST}" = "true" ]; then
    log_section "${LOG_SECTION_PREFIX} Running pipeline self-test"

    TEST_DATA_DIR="${PIPELINE_DIR}/data"
    if [ ! -d "${TEST_DATA_DIR}" ]; then
        log_warn "  Test data directory not found: ${TEST_DATA_DIR} — skipping self-test"
    else
        TEST_OUT_DIR="${OUTPUT_DIR}/test_run"
        mkdir -p "${TEST_OUT_DIR}"

        log_info "  Test input   : ${TEST_DATA_DIR}"
        log_info "  Test output  : ${TEST_OUT_DIR}"

        # Find test input files
        TEST_VCF=$(find "${TEST_DATA_DIR}" -name "*.vcf.gz" | head -1)
        TEST_PHENO=$(find "${TEST_DATA_DIR}" -name "*.tsv" | grep -i pheno | head -1 || \
                     find "${TEST_DATA_DIR}" -name "*.tsv" | head -1)
        TEST_COV=$(find "${TEST_DATA_DIR}" -name "*.tsv" | grep -i cov | head -1 || \
                   find "${TEST_DATA_DIR}" -name "*.tsv" | tail -1)

        if [ -z "${TEST_VCF}" ] || [ -z "${TEST_PHENO}" ] || [ -z "${TEST_COV}" ]; then
            log_warn "  Could not auto-detect test input files — skipping self-test"
            log_warn "  (Expected .vcf.gz + phenotype .tsv + covariate .tsv in ${TEST_DATA_DIR})"
        else
            log_info "  Test VCF     : ${TEST_VCF}"
            log_info "  Test phenotype: ${TEST_PHENO}"
            log_info "  Test covariate: ${TEST_COV}"

            TEST_CMD="${NXF_VER:+NXF_VER=${NXF_VER}} nextflow run ${PIPELINE_DIR}/mvgwas.nf"
            TEST_CMD+=" --geno ${TEST_VCF}"
            TEST_CMD+=" --pheno ${TEST_PHENO}"
            TEST_CMD+=" --cov ${TEST_COV}"
            TEST_CMD+=" --dir ${TEST_OUT_DIR}"
            TEST_CMD+=" --out test_result.tsv"
            TEST_CMD+=" --l 500"
            TEST_CMD+=" ${NXF_CONTAINER_FLAG}"
            TEST_CMD+=" -work-dir ${TEST_OUT_DIR}/.nf_work"

            log_info "  Running: ${TEST_CMD}"

            if [ "${DRY_RUN}" = "true" ]; then
                log_info "  [DRY-RUN] Would run test"
            else
                cd "${TEST_OUT_DIR}"
                eval "${TEST_CMD}" 2>&1 | tee -a "${PIPELINE_LOG}" \
                || die "Pipeline self-test FAILED. Check logs for details."

                # Verify output
                TEST_RESULT="${TEST_OUT_DIR}/result/test_result.tsv"
                if [ -f "${TEST_RESULT}" ] && [ -s "${TEST_RESULT}" ]; then
                    N_LINES=$(wc -l < "${TEST_RESULT}")
                    log_ok "  Self-test PASSED — result: ${TEST_RESULT} (${N_LINES} lines)"
                    record_stat "step1_selftest" "PASSED"
                else
                    die "Self-test output not found or empty: ${TEST_RESULT}"
                fi
            fi
        fi
    fi
else
    log_info "  Self-test skipped (--skip-test / --no-test)"
    record_stat "step1_selftest" "SKIPPED"
fi

log_ok "Step I complete — environment is ready."
