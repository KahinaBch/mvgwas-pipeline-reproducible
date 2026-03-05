#!/usr/bin/env python3
# =============================================================================
# scripts/02_sex_stratify.py  —  Detect sex coding and split sample ID lists
# =============================================================================
"""
Reads the (already-filtered) covariate TSV, auto-detects the sex/gender column
and its coding, then writes separate sample-ID lists for males and females.

Supported codings (auto-detected):
  1/2    (PLINK standard: 1=male, 2=female)
  M/F    (or m/f, any case)
  male/female  (any case)
  0/1    (non-standard; warns user; assumes 1=male, 0=female)

Override with --male-code / --female-code if auto-detection fails.
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser(description="Split samples by sex from covariate file")
parser.add_argument("--cov-in",        required=True, help="Filtered covariate TSV (or .gz)")
parser.add_argument("--sex-col",       default=None,  help="Column name for sex (auto-detected if omitted)")
parser.add_argument("--male-code",     default=None,  help="Override male code value")
parser.add_argument("--female-code",   default=None,  help="Override female code value")
parser.add_argument("--samples-male",   required=True, help="Output: male sample IDs (one per line)")
parser.add_argument("--samples-female", required=True, help="Output: female sample IDs (one per line)")
parser.add_argument("--report",        required=True, help="Output: text report path")
args = parser.parse_args()


def open_tsv(path):
    if str(path).endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path)


# ── Load covariate file ────────────────────────────────────────────────────────
with open_tsv(args.cov_in) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fieldnames = reader.fieldnames or []
    rows = list(reader)

if not rows:
    print("ERROR: Covariate file is empty.", file=sys.stderr)
    sys.exit(1)

# ── Auto-detect sex column ────────────────────────────────────────────────────
sex_col = args.sex_col
if not sex_col:
    # Exact-match candidates first, then partial
    for candidate in ["sex", "Sex", "SEX", "gender", "Gender", "GENDER"]:
        if candidate in fieldnames:
            sex_col = candidate
            break
    if not sex_col:
        for f in fieldnames:
            if "sex" in f.lower() or "gender" in f.lower():
                sex_col = f
                break

if not sex_col:
    print(
        f"ERROR: Could not find a sex/gender column. "
        f"Available columns: {', '.join(fieldnames)}. "
        f"Use --sex-col to specify it explicitly.",
        file=sys.stderr,
    )
    sys.exit(1)

print(f"[INFO] Sex column identified: '{sex_col}'")

# ── Get unique values ─────────────────────────────────────────────────────────
vals = [r[sex_col] for r in rows if r[sex_col] not in ("", "NA", "na", "NaN", "nan", ".")]
val_counts = Counter(vals)
unique_vals = sorted(val_counts.keys())
print(f"[INFO] Unique values in '{sex_col}': {dict(val_counts)}")

if len(unique_vals) < 2:
    print(
        f"ERROR: Expected at least 2 distinct sex values, found only: {unique_vals}",
        file=sys.stderr,
    )
    sys.exit(1)

# ── Auto-detect coding ────────────────────────────────────────────────────────
male_code   = args.male_code
female_code = args.female_code

if not male_code or not female_code:
    uv_lower = {v.lower(): v for v in unique_vals}

    if "1" in unique_vals and "2" in unique_vals:
        male_code, female_code = "1", "2"
        print("[INFO] Detected coding: 1=male, 2=female (PLINK standard)")

    elif "m" in uv_lower and "f" in uv_lower:
        male_code   = uv_lower["m"]
        female_code = uv_lower["f"]
        print(f"[INFO] Detected coding: {male_code}=male, {female_code}=female")

    elif "male" in uv_lower and "female" in uv_lower:
        # "female" contains "male" so check female first
        male_code   = uv_lower["male"]
        female_code = uv_lower["female"]
        print(f"[INFO] Detected coding: '{male_code}'=male, '{female_code}'=female")

    elif "0" in unique_vals and "1" in unique_vals:
        male_code, female_code = "1", "0"
        print(
            "[WARN] Ambiguous 0/1 coding — assuming 1=male, 0=female. "
            "Override with SEX_MALE_CODE and SEX_FEMALE_CODE in config if wrong.",
            file=sys.stderr,
        )

    else:
        print(
            f"ERROR: Cannot auto-detect sex coding from values: {unique_vals}. "
            f"Set SEX_MALE_CODE and SEX_FEMALE_CODE in your config file.",
            file=sys.stderr,
        )
        sys.exit(1)

# ── Split sample IDs ──────────────────────────────────────────────────────────
id_col = fieldnames[0]   # first column = sample ID

male_ids   = [r[id_col] for r in rows if r[sex_col] == male_code]
female_ids = [r[id_col] for r in rows if r[sex_col] == female_code]
other_ids  = [r[id_col] for r in rows if r[sex_col] not in (male_code, female_code)]

print(f"[INFO] Males: {len(male_ids)}, Females: {len(female_ids)}, Other/unknown: {len(other_ids)}")

if len(male_ids) == 0:
    print(f"ERROR: No male samples found (code='{male_code}'). Check SEX_MALE_CODE.", file=sys.stderr)
    sys.exit(1)
if len(female_ids) == 0:
    print(f"ERROR: No female samples found (code='{female_code}'). Check SEX_FEMALE_CODE.", file=sys.stderr)
    sys.exit(1)

# ── Write outputs ─────────────────────────────────────────────────────────────
Path(args.samples_male).write_text("\n".join(male_ids) + "\n")
Path(args.samples_female).write_text("\n".join(female_ids) + "\n")

report = [
    "=== Sex Stratification Report ===",
    f"Covariate file  : {args.cov_in}",
    f"Sex column      : {sex_col}",
    f"Male code       : {male_code}",
    f"Female code     : {female_code}",
    f"",
    f"Male samples    : {len(male_ids)}",
    f"Female samples  : {len(female_ids)}",
    f"Other/unknown   : {len(other_ids)}",
    "",
    "Value counts:",
    *[f"  {v}: {c}" for v, c in sorted(val_counts.items())],
    "",
    "Male sample IDs (first 5):",
    *[f"  {i}" for i in male_ids[:5]],
    "Female sample IDs (first 5):",
    *[f"  {i}" for i in female_ids[:5]],
]
Path(args.report).write_text("\n".join(report) + "\n")
print(f"[OK] Sex stratification complete.")
