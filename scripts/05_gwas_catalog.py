#!/usr/bin/env python3
"""
05_gwas_catalog.py  —  GWAS Catalog mapping for mvGWAS top SNPs

For each stratum's rsID-annotated top SNPs:
  1. Cluster SNPs into loci (same CHR, within --clump-kb of each other)
  2. Pick lead SNP per locus (lowest p-value)
  3. Query EBI GWAS Catalog REST API for exact SNP → direct associations
  4. Query gene-region context (±--window-kb around lead SNP) → nearby catalog hits
  5. Write a Markdown report + TSV summary to --output-dir

Usage:
  python3 05_gwas_catalog.py \\
      --top-snps   results/top_1000_snps_rsid.tsv \\
      --output-dir gwas_catalog/combined \\
      --run-name   my_run \\
      --stratum    combined \\
      --pval-col   P_manta \\
      --top-n      20 \\
      --clump-kb   500 \\
      --window-kb  500
"""

import argparse
import csv
import json
import math
import os
import sys
import time
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

# EBI GWAS Catalog REST endpoints
_CATALOG_SNP_ASSOC = (
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms"
    "/{rsid}/associations?projection=associationBySnp&size=50&page=0"
)
_CATALOG_REGION = (
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms"
    "/search/findByChromosomeBetweenLocations"
    "?chrom={chrom}&startLocation={start}&endLocation={end}&size=200&page=0"
)
API_PAUSE = 0.6   # seconds between requests (be polite to EBI servers)


# ──────────────────────────────────────────────────────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top-snps",   required=True,
                    help="rsID-annotated top SNPs TSV (output of Step IV)")
    ap.add_argument("--output-dir", required=True,
                    help="Directory to write the markdown report and TSV")
    ap.add_argument("--run-name",   required=True)
    ap.add_argument("--stratum",    default="combined")
    ap.add_argument("--pval-col",   default="P_manta",
                    help="Column name for the p-value to sort by")
    ap.add_argument("--top-n",      type=int, default=20,
                    help="Number of top lead loci to query (default: 20)")
    ap.add_argument("--clump-kb",   type=float, default=500,
                    help="Locus clumping window in kb (default: 500)")
    ap.add_argument("--window-kb",  type=float, default=500,
                    help="Region context window in kb around each lead SNP (default: 500)")
    return ap.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
def api_get(url: str) -> dict | None:
    """GET JSON from url with retries; returns None on HTTP 404 or network error."""
    time.sleep(API_PAUSE)
    req = Request(url, headers={"Accept": "application/json",
                                "User-Agent": "mvgwas-pipeline/1.0 (research)"})
    for attempt in range(3):
        try:
            with urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except HTTPError as e:
            if e.code == 404:
                return None
            if e.code in (429, 503):
                wait = (attempt + 1) * 10
                print(f"    Rate-limited (HTTP {e.code}), waiting {wait}s …",
                      file=sys.stderr, flush=True)
                time.sleep(wait)
                continue
            print(f"    HTTP {e.code} for {url}", file=sys.stderr)
            return None
        except URLError as e:
            print(f"    URL error: {e}", file=sys.stderr)
            if attempt < 2:
                time.sleep(5)
                continue
            return None
    return None


def check_connectivity() -> bool:
    """Quick check that the EBI GWAS Catalog API is reachable."""
    try:
        urlopen("https://www.ebi.ac.uk/gwas/rest/api", timeout=10)
        return True
    except Exception:
        return False


# ──────────────────────────────────────────────────────────────────────────────
def _detect_col(header: list[str], candidates: list[str]) -> str | None:
    upper = {h.upper(): h for h in header}
    for c in candidates:
        if c.upper() in upper:
            return upper[c.upper()]
    return None


def load_top_snps(path: str, pval_col: str) -> tuple[list, str | None, str | None, str | None]:
    """Read TSV, return rows sorted ascending by pval_col; also return detected column names."""
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        chr_col = _detect_col(header, ["CHR", "CHROM", "#CHROM", "chromosome"])
        pos_col = _detect_col(header, ["POS", "BP", "POSITION", "start"])
        id_col  = _detect_col(header, ["ID", "RSID", "RS_ID", "SNP", "variant_id"])
        # Detect the actual pval column name (case-insensitive partial match as fallback)
        pval_col_actual = _detect_col(header, [pval_col]) or \
            next((h for h in header if pval_col.lower() in h.lower()), None)
        if pval_col_actual is None:
            print(f"ERROR: p-value column '{pval_col}' not found in {path}.\n"
                  f"Available columns: {header}", file=sys.stderr)
            sys.exit(1)
        for row in reader:
            try:
                p = float(row.get(pval_col_actual, "") or "nan")
            except ValueError:
                p = float("nan")
            row["_pval"] = p
            rows.append(row)
    rows = [r for r in rows if not math.isnan(r["_pval"]) and r["_pval"] > 0]
    rows.sort(key=lambda r: r["_pval"])
    return rows, chr_col, pos_col, id_col


def clump_loci(rows: list, chr_col: str, pos_col: str,
               clump_kb: float) -> list[dict]:
    """Greedy clumping: group rows into loci (same chr, within clump_kb bp)."""
    loci: list[dict] = []
    for row in rows:
        chrom = (row.get(chr_col, "") or "").replace("chr", "").strip()
        try:
            pos = int(row.get(pos_col, 0) or 0)
        except (ValueError, TypeError):
            continue
        placed = False
        for locus in loci:
            if locus["chr"] == chrom and abs(pos - locus["pos_max"]) <= clump_kb * 1000:
                locus["snps"].append(row)
                locus["pos_max"] = max(locus["pos_max"], pos)
                placed = True
                break
        if not placed:
            loci.append({"chr": chrom, "pos_max": pos, "snps": [row]})
    for locus in loci:
        locus["snps"].sort(key=lambda r: r["_pval"])
        locus["lead"] = locus["snps"][0]
    loci.sort(key=lambda l: l["lead"]["_pval"])
    return loci


# ──────────────────────────────────────────────────────────────────────────────
def query_direct(rsid: str) -> list[dict]:
    """Query GWAS Catalog for exact SNP associations. Returns list of assoc dicts."""
    if not rsid or not rsid.startswith("rs"):
        return []
    data = api_get(_CATALOG_SNP_ASSOC.format(rsid=rsid))
    if data is None:
        return []
    try:
        assocs = data["_embedded"]["associations"]
    except (KeyError, TypeError):
        return []
    results = []
    for a in assocs:
        try:
            efo = a.get("efoTraits") or a.get("reportedGenes") or []
            trait = efo[0].get("trait", "?") if efo else "?"
        except (IndexError, TypeError, KeyError):
            trait = "?"
        pval  = a.get("pvalue", "?")
        study = a.get("study") or {}
        pub   = study.get("publicationInfo") or {}
        pmid  = pub.get("pubmedId", "?")
        gcst  = study.get("accessionId", "?")
        results.append({"trait": trait, "pvalue": pval, "pmid": pmid, "gcst": gcst})
    return results


def query_region(chrom: str, pos: int, window_kb: float) -> list[dict]:
    """Query GWAS Catalog for any variants in ±window_kb window."""
    start = max(1, pos - int(window_kb * 1000))
    end   = pos + int(window_kb * 1000)
    data  = api_get(_CATALOG_REGION.format(chrom=chrom, start=start, end=end))
    if data is None:
        return []
    try:
        snps_raw = data["_embedded"]["singleNucleotidePolymorphisms"]
    except (KeyError, TypeError):
        return []
    hits = []
    for snp in snps_raw:
        rs    = snp.get("rsId", "?")
        locs  = snp.get("locations") or []
        bp    = locs[0].get("chromosomePosition", "?") if locs else "?"
        n_ass = len(snp.get("associations", []))
        hits.append({"rsid": rs, "pos": bp, "n_assoc": n_ass})
    return hits


def _fmt_p(pval) -> str:
    try:
        v = float(pval)
        if v <= 0:
            return "?"
        exp = int(math.floor(math.log10(v)))
        m   = v / 10 ** exp
        return f"{m:.2f} × 10<sup>{exp}</sup>"
    except Exception:
        return str(pval)


# ──────────────────────────────────────────────────────────────────────────────
def run(args):
    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.isfile(args.top_snps):
        print(f"ERROR: top SNPs file not found: {args.top_snps}", file=sys.stderr)
        sys.exit(1)

    # ── Check connectivity ────────────────────────────────────────────────────
    print("Checking EBI GWAS Catalog connectivity …", flush=True)
    if not check_connectivity():
        print("WARNING: EBI GWAS Catalog appears unreachable. "
              "Writing a placeholder report.", file=sys.stderr)
        _write_offline_report(args)
        return

    # ── Load SNPs ─────────────────────────────────────────────────────────────
    print(f"Loading top SNPs from: {args.top_snps}", flush=True)
    rows, chr_col, pos_col, id_col = load_top_snps(args.top_snps, args.pval_col)
    if not rows:
        print("ERROR: no valid rows found (check --pval-col)", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(rows)} rows loaded  |  CHR={chr_col}  POS={pos_col}  ID={id_col}",
          flush=True)

    # ── Clump + take top-N loci ───────────────────────────────────────────────
    # Use extra headroom (5×top_n) before clumping so we get top_n distinct loci
    loci = clump_loci(rows[:min(args.top_n * 10, len(rows))],
                      chr_col or "CHR", pos_col or "POS", args.clump_kb)
    loci = loci[:args.top_n]
    print(f"  {len(loci)} distinct loci (clumped at {args.clump_kb:.0f} kb)", flush=True)

    # ── Query GWAS Catalog ────────────────────────────────────────────────────
    results = []
    for i, locus in enumerate(loci):
        lead  = locus["lead"]
        rsid  = (lead.get(id_col, "") if id_col else "").strip()
        chrom = ((lead.get(chr_col, "") if chr_col else "")
                 .replace("chr", "").strip())
        try:
            pos = int(lead.get(pos_col, 0) or 0)
        except (ValueError, TypeError):
            pos = 0
        pval = lead["_pval"]

        print(f"  [{i+1:2d}/{len(loci)}] {rsid:<15}  CHR{chrom}:{pos:,}  "
              f"P={pval:.2e}", flush=True)

        direct      = query_direct(rsid)
        region_hits = query_region(chrom, pos, args.window_kb)

        results.append({
            "rsid":         rsid,
            "chr":          chrom,
            "pos":          pos,
            "pval":         pval,
            "n_locus":      len(locus["snps"]),
            "direct":       direct,
            "region_hits":  region_hits,
        })

    # ── Write TSV ─────────────────────────────────────────────────────────────
    tsv_path = os.path.join(args.output_dir,
                            f"gwas_catalog_{args.stratum}.tsv")
    with open(tsv_path, "w") as fh:
        fh.write("rsID\tCHR\tPOS\tP_GWAS\tN_locus_SNPs\t"
                 "N_direct_hits\tN_region_SNPs\tTop_trait\tTop_PMID\tTop_GCST\n")
        for r in results:
            d = r["direct"]
            top_trait = d[0]["trait"] if d else "Novel"
            top_pmid  = d[0]["pmid"]  if d else "—"
            top_gcst  = d[0]["gcst"]  if d else "—"
            fh.write(
                f"{r['rsid']}\t{r['chr']}\t{r['pos']}\t{r['pval']:.4e}\t"
                f"{r['n_locus']}\t{len(d)}\t{len(r['region_hits'])}\t"
                f"{top_trait}\t{top_pmid}\t{top_gcst}\n"
            )
    print(f"  TSV summary → {tsv_path}", flush=True)

    # ── Write Markdown report ─────────────────────────────────────────────────
    _write_md_report(args, results)


def _write_offline_report(args):
    """Minimal placeholder when the API is not reachable."""
    md_path = os.path.join(args.output_dir,
                           f"gwas_catalog_{args.stratum}.md")
    with open(md_path, "w") as fh:
        fh.write(f"# GWAS Catalog Mapping — {args.run_name} ({args.stratum})\n\n")
        fh.write(f"> **WARNING:** EBI GWAS Catalog API was not reachable on "
                 f"{time.strftime('%Y-%m-%d')}. Re-run with internet access.\n\n")
    print(f"  Placeholder report → {md_path}", flush=True)


def _write_md_report(args, results: list[dict]):
    """Write full Markdown report."""
    md_path = os.path.join(args.output_dir,
                           f"gwas_catalog_{args.stratum}.md")
    n_novel   = sum(1 for r in results if not r["direct"])
    n_catalog = len(results) - n_novel

    with open(md_path, "w") as fh:
        fh.write(f"# GWAS Catalog Mapping — {args.run_name} ({args.stratum})\n\n")
        fh.write(f"**Date:** {time.strftime('%Y-%m-%d')}  \n")
        fh.write(f"**Stratum:** {args.stratum}  \n")
        fh.write(f"**Input file:** `{os.path.basename(args.top_snps)}`  \n")
        fh.write(f"**Lead SNPs queried:** {len(results)} "
                 f"(top {args.top_n} distinct loci, "
                 f"clumped at {args.clump_kb:.0f} kb)  \n\n")
        fh.write("---\n\n")

        # ── Method ────────────────────────────────────────────────────────────
        fh.write("## Method\n\n")
        fh.write(
            f"Top {args.top_n} lead SNPs (one per locus, clumped at "
            f"{args.clump_kb:.0f} kb) were queried against the "
            f"[EBI GWAS Catalog](https://www.ebi.ac.uk/gwas) REST API:\n\n"
        )
        fh.write(
            "1. **Direct SNP lookup** via "
            "`/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations`  \n"
        )
        fh.write(
            f"2. **Region context** via the region endpoint (±{args.window_kb:.0f} kb "
            "around each lead SNP) — captures prior associations at nearby loci.\n\n"
        )
        fh.write("---\n\n")

        # ── Lead SNPs table ────────────────────────────────────────────────────
        fh.write("## Lead SNPs\n\n")
        fh.write(f"| # | rsID | CHR | POS | {args.pval_col} | "
                 f"Locus SNPs | Direct hits | ±{args.window_kb:.0f} kb catalog SNPs |\n")
        fh.write("|---|------|-----|-----|------------|-----------|------------|---|\n")
        for i, r in enumerate(results, 1):
            d_str = f"✅ {len(r['direct'])}" if r["direct"] else "❌ None"
            fh.write(
                f"| {i} | {r['rsid']} | {r['chr']} | {r['pos']:,} | "
                f"{r['pval']:.2e} | {r['n_locus']} | {d_str} | "
                f"{len(r['region_hits'])} |\n"
            )
        fh.write("\n---\n\n")

        # ── Results by locus ──────────────────────────────────────────────────
        fh.write("## Results by Locus\n\n")
        fh.write(
            f"**{n_novel}/{len(results)} lead SNPs are novel** "
            "(no direct GWAS Catalog entry).  \n"
        )
        if n_catalog:
            fh.write(
                f"**{n_catalog}/{len(results)} lead SNPs have prior catalog associations.**\n\n"
            )
        else:
            fh.write("\n")

        for i, r in enumerate(results, 1):
            fh.write(
                f"### Locus {i}: {r['rsid']}  —  "
                f"CHR{r['chr']}:{r['pos']:,}\n\n"
            )
            fh.write(f"- **{args.pval_col}**: `{r['pval']:.4e}`  \n")
            fh.write(f"- **Locus size**: {r['n_locus']} SNP(s) within "
                     f"{args.clump_kb:.0f} kb  \n")
            fh.write(
                f"- **Nearby GWAS Catalog SNPs** "
                f"(±{args.window_kb:.0f} kb): {len(r['region_hits'])}  \n\n"
            )
            if r["direct"]:
                fh.write(
                    f"**Direct associations ({len(r['direct'])}):**\n\n"
                )
                fh.write("| Trait | P-value | PMID | Study |\n")
                fh.write("|-------|---------|------|-------|\n")
                for a in r["direct"][:25]:
                    fh.write(
                        f"| {a['trait']} | {_fmt_p(a['pvalue'])} | "
                        f"{a['pmid']} | {a['gcst']} |\n"
                    )
                fh.write("\n")
            else:
                fh.write(
                    "**Direct associations: ❌ Novel** — this exact SNP "
                    "has no prior GWAS Catalog entry.\n\n"
                )
            if r["region_hits"]:
                top_region = sorted(
                    r["region_hits"],
                    key=lambda x: x.get("n_assoc", 0),
                    reverse=True
                )[:5]
                fh.write(
                    f"**Top nearby catalog variants (±{args.window_kb:.0f} kb):**  \n"
                )
                for h in top_region:
                    fh.write(f"- {h['rsid']} (pos {h['pos']}, "
                             f"{h['n_assoc']} assoc)  \n")
                fh.write("\n")

        # ── Summary table ─────────────────────────────────────────────────────
        fh.write("---\n\n")
        fh.write("## Summary Table\n\n")
        fh.write(
            f"| Locus | rsID | CHR:POS | {args.pval_col} | "
            "Direct hits | Top trait |\n"
        )
        fh.write("|-------|------|---------|------------|-------------|----------|\n")
        for r in results:
            top_trait = r["direct"][0]["trait"] if r["direct"] else "Novel"
            fh.write(
                f"| chr{r['chr']}:{r['pos']//1_000_000:.0f} Mb | "
                f"{r['rsid']} | {r['chr']}:{r['pos']:,} | "
                f"{r['pval']:.2e} | {len(r['direct'])} | {top_trait} |\n"
            )
        fh.write("\n")
        fh.write(
            f"*Queried: EBI GWAS Catalog (https://www.ebi.ac.uk/gwas), "
            f"accessed {time.strftime('%Y-%m-%d')}*  \n"
        )
        fh.write(
            f"*Genome assembly: GRCh38, "
            f"clumping window: {args.clump_kb:.0f} kb, "
            f"region window: {args.window_kb:.0f} kb*\n"
        )

    print(f"  Markdown report → {md_path}", flush=True)
    print("Done.")


if __name__ == "__main__":
    run(parse_args())
