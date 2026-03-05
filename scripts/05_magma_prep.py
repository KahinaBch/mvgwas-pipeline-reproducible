#!/usr/bin/env python3
"""
05_magma_prep.py  —  Prepare MAGMA input files from a merged mvGWAS results file.

Produces:
  --snp-loc  : SNP  CHR  BP        (for magma --annotate)
  --pval-out : SNP  P               (for magma --pval)

Only rows with a proper rsID (starts with 'rs') and a valid numeric p-value
are written.  Chromosome prefixes ('chr') are stripped automatically.

Usage:
  python3 05_magma_prep.py \\
      --merged-results results/mvgwas_merged.tsv \\
      --pval-col       P_manta \\
      --snp-loc        magma/input/snp_loc.txt \\
      --pval-out       magma/input/pval.txt
"""

import argparse
import csv
import os
import sys


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-results", required=True,
                    help="Merged GWAS TSV (must have rsID, CHR, POS, and p-value columns)")
    ap.add_argument("--pval-col",       default="P_manta",
                    help="Name of the p-value column (default: P_manta)")
    ap.add_argument("--snp-loc",        required=True,
                    help="Output path for the SNP-location file (SNP CHR BP)")
    ap.add_argument("--pval-out",       required=True,
                    help="Output path for the p-value file (SNP P)")
    return ap.parse_args()


def _detect(header: list[str], candidates: list[str]) -> str | None:
    up = {h.upper(): h for h in header}
    for c in candidates:
        if c.upper() in up:
            return up[c.upper()]
    return None


def run(args):
    os.makedirs(os.path.dirname(os.path.abspath(args.snp_loc)),  exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.pval_out)), exist_ok=True)

    written = skip_rsid = skip_pval = 0

    with open(args.merged_results) as fin, \
         open(args.snp_loc, "w") as fsnp, \
         open(args.pval_out, "w") as fpval:

        reader = csv.DictReader(fin, delimiter="\t")
        header = reader.fieldnames or []

        # Detect column names (case-insensitive, with fallbacks)
        chr_col  = _detect(header, ["CHR", "CHROM", "#CHROM", "chromosome"])
        pos_col  = _detect(header, ["POS", "BP", "POSITION", "start"])
        id_col   = _detect(header, ["ID", "RSID", "RS_ID", "SNP", "variant_id"])
        pval_col = _detect(header, [args.pval_col]) or \
            next((h for h in header if args.pval_col.lower() in h.lower()), None)

        if pval_col is None:
            print(
                f"ERROR: p-value column '{args.pval_col}' not found.\n"
                f"Available columns: {header}",
                file=sys.stderr,
            )
            sys.exit(1)

        print(f"Column map:  CHR={chr_col}  POS={pos_col}  "
              f"ID={id_col}  P={pval_col}", flush=True)

        # Write headers
        fsnp.write("SNP\tCHR\tBP\n")
        fpval.write("SNP\tP\n")

        for row in reader:
            rsid = (row.get(id_col, "") if id_col else "").strip()
            if not rsid or not rsid.lower().startswith("rs"):
                skip_rsid += 1
                continue

            praw = (row.get(pval_col, "") if pval_col else "").strip()
            if not praw or praw in (".", "NA", "nan", "NaN", ""):
                skip_pval += 1
                continue
            try:
                float(praw)
            except ValueError:
                skip_pval += 1
                continue

            chrom = (row.get(chr_col, "") if chr_col else "").replace("chr", "").strip()
            bp    = (row.get(pos_col, "") if pos_col else "").strip()

            fsnp.write(f"{rsid}\t{chrom}\t{bp}\n")
            fpval.write(f"{rsid}\t{praw}\n")
            written += 1

            if written % 500_000 == 0:
                print(f"  {written:,} SNPs written …", flush=True)

    print(f"\nDone.")
    print(f"  Written  (rsID + p-value) : {written:,}")
    print(f"  Skipped  (no/invalid rsID): {skip_rsid:,}")
    print(f"  Skipped  (no p-value)     : {skip_pval:,}")
    print(f"  SNP-loc  → {args.snp_loc}")
    print(f"  P-values → {args.pval_out}")


if __name__ == "__main__":
    run(parse_args())
