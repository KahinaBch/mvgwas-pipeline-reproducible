#!/usr/bin/env python3
"""
scripts/02_filter_individuals.py
=================================
Filter phenotype and covariate TSV files to a set of sample IDs,
preserving the original column order and header.

Called by 02_check_inputs.sh.
"""

import argparse
import sys
import os
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Filter phenotype and covariate files to a sample ID intersection."
    )
    p.add_argument("--pheno-in",  required=True,  help="Input phenotype TSV (or .tsv.gz)")
    p.add_argument("--cov-in",    required=True,  help="Input covariate TSV (or .tsv.gz)")
    p.add_argument("--samples",   required=True,  help="File with one sample ID per line")
    p.add_argument("--pheno-out", required=True,  help="Output phenotype TSV")
    p.add_argument("--cov-out",   required=True,  help="Output covariate TSV")
    p.add_argument("--report",    default=None,   help="Write a text summary report here")
    p.add_argument("--verbose",   action="store_true")
    return p.parse_args()


def vprint(msg, verbose):
    if verbose:
        print(f"  [filter] {msg}", flush=True)


def read_tsv(path, verbose):
    """Read a TSV or TSV.GZ into a list of lines (strings)."""
    import gzip
    vprint(f"Reading {path}", verbose)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.readlines()
    return lines


def filter_tsv(in_path, keep_set, out_path, verbose):
    """
    Filter a TSV to rows whose first-column value is in keep_set.
    Returns (n_total, n_kept, n_dropped, list_of_missing_ids).
    """
    lines = read_tsv(in_path, verbose)
    if not lines:
        raise ValueError(f"Empty file: {in_path}")

    header = lines[0]
    data   = lines[1:]
    n_total = len(data)

    kept   = []
    found  = set()
    for line in data:
        sid = line.split("\t")[0].strip()
        if sid in keep_set:
            kept.append(line)
            found.add(sid)

    n_kept = len(kept)
    n_dropped = n_total - n_kept

    # IDs in keep_set but not in this file
    missing = sorted(keep_set - found)

    vprint(f"  -> {n_total} rows in, {n_kept} kept, {n_dropped} dropped", verbose)
    if missing:
        vprint(f"  -> {len(missing)} IDs from intersection not found in file", verbose)

    with open(out_path, "w") as fh:
        fh.write(header)
        for line in kept:
            fh.write(line)

    return n_total, n_kept, n_dropped, missing


def main():
    args = parse_args()

    # Load sample intersection
    with open(args.samples, "r") as fh:
        keep_set = set(line.strip() for line in fh if line.strip())

    n_keep = len(keep_set)
    vprint(f"Target sample set: {n_keep} IDs", args.verbose)

    # Filter phenotype
    ph_total, ph_kept, ph_dropped, ph_missing = filter_tsv(
        args.pheno_in, keep_set, args.pheno_out, args.verbose
    )

    # Filter covariate
    cv_total, cv_kept, cv_dropped, cv_missing = filter_tsv(
        args.cov_in, keep_set, args.cov_out, args.verbose
    )

    # Consistency check: both filtered files must have the same number of rows
    if ph_kept != cv_kept:
        print(
            f"ERROR: phenotype filtered={ph_kept} rows but covariate filtered={cv_kept} rows. "
            "This should not happen if the intersection was computed correctly.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Report
    lines = [
        "=" * 60,
        "  Individual filtering report",
        f"  Target IDs  : {n_keep}",
        "",
        f"  Phenotype file: {args.pheno_in}",
        f"    Total rows  : {ph_total}",
        f"    Kept        : {ph_kept}",
        f"    Dropped     : {ph_dropped}",
    ]
    if ph_missing:
        lines.append(f"    Not in file : {len(ph_missing)} IDs (first 5: {ph_missing[:5]})")

    lines += [
        "",
        f"  Covariate file: {args.cov_in}",
        f"    Total rows  : {cv_total}",
        f"    Kept        : {cv_kept}",
        f"    Dropped     : {cv_dropped}",
    ]
    if cv_missing:
        lines.append(f"    Not in file : {len(cv_missing)} IDs (first 5: {cv_missing[:5]})")

    lines += [
        "",
        f"  Final sample size : {ph_kept}",
        "=" * 60,
    ]

    report_text = "\n".join(lines)
    print(report_text)

    if args.report:
        with open(args.report, "w") as fh:
            fh.write(report_text + "\n")
        vprint(f"Report written: {args.report}", args.verbose)

    if ph_kept == 0:
        print("ERROR: No samples retained after filtering!", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
