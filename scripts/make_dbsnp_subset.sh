#!/usr/bin/env bash
# =============================================================================
# scripts/make_dbsnp_subset.sh  —  small local dbSNP VCF for Step IV (rsIDs)
# =============================================================================
# Step IV annotates the top-N SNPs of each stratum with rsIDs from dbSNP. The
# full NCBI dbSNP VCF is ~30 GB, and the remote fallback needs internet access
# on the node running Step IV (HPC compute nodes often have none) and is very
# slow (one HTTPS request per position). This script, run once on a machine
# WITH internet access (e.g. a login node) after Step III, extracts only the
# positions Step IV will look up and writes a small, indexed VCF with contigs
# renamed to 1..22 like the analysis VCF. Point DBSNP_VCF in the config at it.
#
# Usage:
#   bash scripts/make_dbsnp_subset.sh OUT.vcf.gz TOP_N merged1.tsv [merged2.tsv ...]
#     TOP_N        same value as TOP_N_SNPS in the config (default there: 1000)
#     mergedX.tsv  <OUTPUT_DIR>/results*/mvgwas_merged.tsv of each stratum/run
#
# Environment:
#   DBSNP     dbSNP VCF (local path or URL). Default: NCBI GRCh38 latest release.
#             GRCh37: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
#             A local copy is much faster than the URL; download it with
#               curl -O <URL> && curl -O <URL>.tbi    (~30 GB, a few minutes)
#   PVAL_COL  p-value column used to rank SNPs (default P_manta, as Step IV)
#
# Example (all three strata of one run):
#   DBSNP=/data/dbsnp/GCF_000001405.40.gz bash scripts/make_dbsnp_subset.sh \
#       /data/dbsnp/dbsnp_top1000.vcf.gz 1000 OUT/results{,_female,_male}/mvgwas_merged.tsv
# =============================================================================
set -euo pipefail
if [ $# -lt 3 ]; then
    sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
OUT="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; TOP_N="$2"; shift 2
DBSNP="${DBSNP:-https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz}"
PVAL_COL="${PVAL_COL:-P_manta}"
for t in bcftools tabix awk sort; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not found in PATH" >&2; exit 1; }
done

MERGED=()
for f in "$@"; do
    [ -s "$f" ] || { echo "ERROR: not found or empty: $f" >&2; exit 1; }
    MERGED+=("$(cd "$(dirname "$f")" && pwd)/$(basename "$f")")
done
case "$DBSNP" in http*|ftp*) ;; *) DBSNP="$(cd "$(dirname "$DBSNP")" && pwd)/$(basename "$DBSNP")" ;; esac

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
cd "$tmp"   # a remote index (.tbi) is downloaded to the working directory

# NCBI contig <-> chromosome (1..22, X, Y); identity if dbSNP already uses 1..22
tabix -l "$DBSNP" | awk '{c=$1; if (c ~ /^NC_0000[0-9][0-9]\./) {k=substr(c,8,2)+0; if(k==23)k="X"; if(k==24)k="Y"} else {k=c; sub(/^chr/,"",k)} print k"\t"c}' > chrmap.tsv
awk -v OFS='\t' '$1!=$2 {print $2, $1}' chrmap.tsv > rename.tsv

# Top-N positions of each merged file, ranked as in 04_postprocess.sh
for f in "${MERGED[@]}"; do
    col=$(head -1 "$f" | tr '\t' '\n' | grep -n "^${PVAL_COL}\$" | cut -d: -f1)
    [ -n "$col" ] || { echo "ERROR: no ${PVAL_COL} column in $f" >&2; exit 1; }
    awk -v col="$col" 'NR>1 && $col!="" && $col!="NA" && $col+0==$col' "$f" \
        | sort -k"${col},${col}"g | awk -v n="$TOP_N" 'NR<=n {c=$1; sub(/^chr/,"",c); print c"\t"$2}'
    echo "  $f: top $TOP_N taken" >&2
done | sort -u > positions.tsv

awk -v OFS='\t' 'NR==FNR{m[$1]=$2; next} ($1 in m){print m[$1], $2, $2}' chrmap.tsv positions.tsv \
    | sort -k1,1 -k2,2n > regions.tsv
echo "Positions to fetch: $(wc -l < regions.tsv)" >&2

if [ -s rename.tsv ]; then
    bcftools view -R regions.tsv -G -Ou "$DBSNP" \
        | bcftools annotate --rename-chrs rename.tsv -x INFO -Ou \
        | bcftools sort -Oz -o "$OUT"
else
    bcftools view -R regions.tsv -G -Ou "$DBSNP" \
        | bcftools annotate -x INFO -Ou \
        | bcftools sort -Oz -o "$OUT"
fi
bcftools index -t -f "$OUT"
echo "Wrote $OUT: $(bcftools index -n "$OUT") dbSNP records" >&2
