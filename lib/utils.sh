#!/usr/bin/env bash
# =============================================================================
# lib/utils.sh  —  Shared utility functions for mvgwas-pipeline
# =============================================================================

# ── Run command with logging and dry-run support ──────────────────────────────
run_cmd() {
    # Usage: run_cmd <description> <cmd> [args...]
    local desc="$1"; shift
    log_info "${desc}"
    if [ "${DRY_RUN:-false}" = "true" ]; then
        log_info "  [DRY-RUN] Would run: $*"
        return 0
    fi
    log_debug "  Executing: $*"
    "$@"
    local rc=$?
    if [ ${rc} -ne 0 ]; then
        die "Command failed (exit ${rc}): $*"
    fi
    return 0
}

# ── Run command, capture output, log it ──────────────────────────────────────
run_cmd_capture() {
    local desc="$1"; shift
    log_info "${desc}"
    if [ "${DRY_RUN:-false}" = "true" ]; then
        log_info "  [DRY-RUN] Would run: $*"
        echo ""
        return 0
    fi
    log_debug "  Executing: $*"
    local out
    out=$("$@" 2>&1) || die "Command failed: $*\nOutput: ${out}"
    echo "${out}"
}

# ── Check that a file exists and is non-empty ─────────────────────────────────
require_file() {
    local path="$1"
    local label="${2:-file}"
    if [ ! -f "${path}" ]; then
        die "${label} not found: ${path}"
    fi
    if [ ! -s "${path}" ]; then
        die "${label} is empty: ${path}"
    fi
    log_ok "${label}: ${path}  ($(wc -l < "${path}") lines, $(du -sh "${path}" | cut -f1))"
}

# ── Check that a file has a .tbi or .csi index ───────────────────────────────
require_index() {
    local path="$1"
    if [ ! -f "${path}.tbi" ] && [ ! -f "${path}.csi" ]; then
        log_warn "No tabix index found for ${path} — attempting to create one..."
        if [[ "${path}" == *.vcf.gz ]]; then
            run_cmd "Indexing ${path}" bcftools index -t "${path}"
        else
            die "Cannot auto-index ${path}: expected .vcf.gz"
        fi
    else
        log_ok "Index OK: ${path}.tbi or ${path}.csi"
    fi
}

# ── Count lines (excluding header) in a TSV/TSV.GZ ───────────────────────────
count_data_rows() {
    local f="$1"
    if [[ "${f}" == *.gz ]]; then
        zcat "${f}" | tail -n +2 | wc -l | tr -d ' '
    else
        tail -n +2 "${f}" | wc -l | tr -d ' '
    fi
}

# ── Get the header of a TSV/TSV.GZ ────────────────────────────────────────────
get_header() {
    local f="$1"
    if [[ "${f}" == *.gz ]]; then
        zcat "${f}" | head -1
    else
        head -1 "${f}"
    fi
}

# ── Pretty-print a number with thousands separator ────────────────────────────
fmt_num() { printf "%'d" "$1"; }

# ── Resolve absolute path (macOS + Linux compatible) ─────────────────────────
abspath() {
    local p="$1"
    if command -v realpath &>/dev/null; then
        realpath "${p}"
    elif command -v greadlink &>/dev/null; then
        greadlink -f "${p}"
    else
        echo "$(cd "$(dirname "${p}")" && pwd)/$(basename "${p}")"
    fi
}

# ── Write a key=value line to the run summary TSV ────────────────────────────
record_stat() {
    local key="$1"
    local val="$2"
    local stats_file="${OUTPUT_DIR}/run_stats.tsv"
    echo -e "${key}\t${val}" >> "${stats_file}"
}

# ── Detect execution environment (slurm / local) ─────────────────────────────
detect_environment() {
    if command -v sbatch &>/dev/null && [ "${USE_SLURM:-false}" = "true" ]; then
        echo "slurm"
    else
        echo "local"
    fi
}

# ── Wait for all background jobs, die if any failed ──────────────────────────
wait_jobs() {
    local desc="${1:-background jobs}"
    log_info "Waiting for ${desc} to finish..."
    local failed=0
    for pid in "${BG_PIDS[@]:-}"; do
        wait "${pid}" || { log_warn "Job PID ${pid} failed"; (( failed++ )) || true; }
    done
    BG_PIDS=()
    if [ "${failed}" -gt 0 ]; then
        die "${failed} ${desc} failed."
    fi
    log_ok "All ${desc} completed successfully."
}
