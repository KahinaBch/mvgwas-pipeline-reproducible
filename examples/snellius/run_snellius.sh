#!/usr/bin/env bash
# =============================================================================
# examples/snellius/run_snellius.sh  —  run mvgwas-pipeline on Snellius (SURF)
# =============================================================================
#   bash   run_snellius.sh check   CONFIG            # login node: steps 1-2
#                                                    # (tools, inputs, sample overlap, strata)
#   sbatch run_snellius.sh run     CONFIG [run.sh options]
#                                                    # steps 1-3 in this job; with USE_SLURM=true
#                                                    # steps 4-6 are queued automatically after
#                                                    # the chromosome + merge jobs
#   bash   run_snellius.sh compare CONFIG MD5FILE    # compare merged results with expected md5s
#
# The modules come from SLURM_SETUP in the config (also used by every job the
# pipeline submits). MD5FILE has one "<md5> <stratum>" line per stratum
# (combined | male | female), e.g. made from an earlier run with
#   md5sum OUT/results/mvgwas_merged.tsv | awk '{print $1, "combined"}'
#
# Snellius notes:
#   - Home is limited to 1M files: put OUTPUT_DIR on project space; Nextflow
#     creates ~10 files per 500-variant chunk.
#   - Compute nodes have no internet: set DBSNP_VCF to a local file (see
#     scripts/make_dbsnp_subset.sh, run on the login node after step 3).
#   - Jobs are billed for at least 1/8 node (24 CPUs) on genoa.
# =============================================================================
#SBATCH --job-name=mvgwas_launch
#SBATCH --output=mvgwas_launch_%j.out
#SBATCH --error=mvgwas_launch_%j.err
#SBATCH --partition=genoa
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail
usage() { sed -n '5,14p' "${SELF:-$0}" | sed 's/^# \{0,1\}//'; exit 2; }

MODE="${1:-}"; CONF="${2:-}"
SELF="${BASH_SOURCE[0]}"
# Submitted with sbatch, this file runs from a spooled copy (no repo next to it):
# then ask SLURM for the real path. Called with bash (also inside another job),
# its own location is used.
if [ ! -f "$(dirname "${SELF}")/../../run.sh" ] && [ -n "${SLURM_JOB_ID:-}" ] \
   && command -v scontrol >/dev/null 2>&1; then
    SELF=$(scontrol show job "${SLURM_JOB_ID}" | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
fi
REPO="${MVGWAS_REPO:-$(cd "$(dirname "${SELF}")/../.." && pwd)}"
[ -f "${REPO}/run.sh" ] || { echo "ERROR: run.sh not found in ${REPO}; set MVGWAS_REPO" >&2; exit 1; }
[ -n "${MODE}" ] && [ -f "${CONF}" ] || usage
CONF="$(cd "$(dirname "${CONF}")" && pwd)/$(basename "${CONF}")"
shift 2

# Environment (modules) from the config's SLURM_SETUP
SETUP=$(bash -c 'set -a; source "$1" >/dev/null 2>&1; printf "%s" "${SLURM_SETUP:-}"' _ "${CONF}")
set +u; eval "${SETUP}"; set -u

compare() {
    local md5file="$1" out status=0
    out=$(bash -c 'source "$1" >/dev/null 2>&1; printf "%s" "${OUTPUT_DIR}"' _ "${CONF}")
    echo "== Merged results vs ${md5file}"
    while read -r want stratum; do
        [ -z "${want}" ] && continue
        local f="${out}/results/mvgwas_merged.tsv"
        [ "${stratum}" != "combined" ] && f="${out}/results_${stratum}/mvgwas_merged.tsv"
        if [ ! -s "${f}" ]; then echo "MISSING    ${stratum}: ${f}"; status=1; continue; fi
        local got; got=$(md5sum "${f}" | cut -d' ' -f1)
        if [ "${got}" = "${want}" ]; then
            echo "IDENTICAL  ${stratum}: ${f}"
        else
            echo "DIFFERENT  ${stratum}: ${f} ($(( $(wc -l < "${f}") - 1 )) variants)"; status=1
        fi
    done < "${md5file}"
    return ${status}
}

cd "${REPO}"
case "${MODE}" in
    check)   bash run.sh --config "${CONF}" --skip-test --steps 1,2 --resume "$@" ;;
    run)     bash run.sh --config "${CONF}" --resume "$@" ;;
    compare) [ -f "${1:-}" ] || usage; compare "$1" ;;
    *)       usage ;;
esac
