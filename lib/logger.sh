#!/usr/bin/env bash
# =============================================================================
# lib/logger.sh  —  Logging utilities for mvgwas-pipeline
# =============================================================================
# Source this file to get coloured, timestamped logging functions.
# All log output goes to STDOUT **and** to ${PIPELINE_LOG} if that variable is
# set (it is set by run.sh before sourcing this file).
# =============================================================================

# ── ANSI colours (disabled when not a TTY) ────────────────────────────────────
if [ -t 1 ]; then
    CLR_RESET='\033[0m'
    CLR_BOLD='\033[1m'
    CLR_CYAN='\033[0;36m'
    CLR_GREEN='\033[0;32m'
    CLR_YELLOW='\033[0;33m'
    CLR_RED='\033[0;31m'
    CLR_MAGENTA='\033[0;35m'
    CLR_BLUE='\033[0;34m'
else
    CLR_RESET=''; CLR_BOLD=''; CLR_CYAN=''; CLR_GREEN=''
    CLR_YELLOW=''; CLR_RED=''; CLR_MAGENTA=''; CLR_BLUE=''
fi

# ── Internal writer ──────────────────────────────────────────────────────────
_write_log() {
    local msg="$1"
    echo -e "${msg}"
    if [ -n "${PIPELINE_LOG:-}" ]; then
        # strip ANSI codes before writing to file
        echo -e "${msg}" | sed 's/\x1b\[[0-9;]*m//g' >> "${PIPELINE_LOG}"
    fi
}

# ── Timestamp helper ─────────────────────────────────────────────────────────
_ts() { date '+%Y-%m-%d %H:%M:%S'; }

# ── Public log functions ──────────────────────────────────────────────────────
log_info()    { _write_log "${CLR_CYAN}[$(_ts)] [INFO ]${CLR_RESET}  $*"; }
log_ok()      { _write_log "${CLR_GREEN}[$(_ts)] [ OK  ]${CLR_RESET}  $*"; }
log_warn()    { _write_log "${CLR_YELLOW}[$(_ts)] [WARN ]${CLR_RESET}  $*"; }
log_error()   { _write_log "${CLR_RED}[$(_ts)] [ERROR]${CLR_RESET}  $*" >&2; }
log_debug()   {
    if [ "${VERBOSE:-false}" = "true" ]; then
        _write_log "${CLR_MAGENTA}[$(_ts)] [DEBUG]${CLR_RESET}  $*"
    fi
}
log_step()    {
    local n="$1"; shift
    _write_log ""
    _write_log "${CLR_BOLD}${CLR_BLUE}[$(_ts)] ════════════════════════════════════════════${CLR_RESET}"
    _write_log "${CLR_BOLD}${CLR_BLUE}[$(_ts)]  STEP ${n}: $*${CLR_RESET}"
    _write_log "${CLR_BOLD}${CLR_BLUE}[$(_ts)] ════════════════════════════════════════════${CLR_RESET}"
    _write_log ""
}
log_section() {
    _write_log ""
    _write_log "${CLR_BOLD}[$(_ts)] ── $* ──────────────────────────────────${CLR_RESET}"
}
log_milestone() {
    _write_log ""
    _write_log "${CLR_GREEN}${CLR_BOLD}[$(_ts)] ✔ MILESTONE: $*${CLR_RESET}"
    _write_log ""
}
log_banner() {
    _write_log ""
    _write_log "${CLR_BOLD}${CLR_CYAN}╔══════════════════════════════════════════════════════╗${CLR_RESET}"
    _write_log "${CLR_BOLD}${CLR_CYAN}║  $*${CLR_RESET}"
    _write_log "${CLR_BOLD}${CLR_CYAN}╚══════════════════════════════════════════════════════╝${CLR_RESET}"
    _write_log ""
}

# ── Checkpoint helpers ────────────────────────────────────────────────────────
# Write a checkpoint file to signal that a step completed successfully.
checkpoint_set() {
    local step="$1"
    local ckpt_dir="${OUTPUT_DIR}/.checkpoints"
    mkdir -p "${ckpt_dir}"
    echo "$(_ts)" > "${ckpt_dir}/step_${step}.done"
    log_debug "Checkpoint written: step_${step}.done"
}

checkpoint_exists() {
    local step="$1"
    [ -f "${OUTPUT_DIR}/.checkpoints/step_${step}.done" ]
}

checkpoint_clear() {
    local step="$1"
    rm -f "${OUTPUT_DIR}/.checkpoints/step_${step}.done"
    log_debug "Checkpoint cleared: step_${step}.done"
}

# ── Elapsed time ──────────────────────────────────────────────────────────────
PIPELINE_START_TS=${PIPELINE_START_TS:-$(date +%s)}

elapsed() {
    local now; now=$(date +%s)
    local diff=$(( now - PIPELINE_START_TS ))
    printf '%02dh %02dm %02ds' $(( diff/3600 )) $(( (diff%3600)/60 )) $(( diff%60 ))
}

# ── die: log error and exit ───────────────────────────────────────────────────
die() {
    log_error "$*"
    log_error "Pipeline aborted after $(elapsed). Check ${PIPELINE_LOG:-stdout} for details."
    exit 1
}

# ── require_cmd: die if command not found ─────────────────────────────────────
require_cmd() {
    local cmd="$1"
    command -v "${cmd}" &>/dev/null || die "Required command not found: ${cmd}"
    log_debug "  command ok: ${cmd} ($(command -v "${cmd}"))"
}
