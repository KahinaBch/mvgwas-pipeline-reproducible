#!/usr/bin/env python3
"""
05_magma_visualize.py  —  MAGMA results visualisation + summary tables

Reads MAGMA output files produced by 05_magma.sh and generates:
  • gene_manhattan.png          — gene-based –log10(P) across chromosomes
  • hallmark_barplot.png        — top 20 Hallmark gene-set enrichments
  • c2_bubbleplot.png           — top 25 C2 curated pathway enrichments
  • c5_barplot.png              — top 20 C5 GO-BP + HPO enrichments
  • top_genes.tsv               — top 50 gene-based associations
  • top_genesets_Hallmark.tsv   — top 50 Hallmark gene sets
  • top_genesets_C2_Curated.tsv — top 50 C2 gene sets
  • top_genesets_C5_GO_HPO.tsv  — top 50 C5 gene sets

Usage:
  python3 05_magma_visualize.py \\
      --output-prefix magma/combined/output/my_run_combined \\
      --plots-dir     magma/combined/plots \\
      --gene-loc      /path/to/NCBI38.gene.loc \\
      --run-name      "my_run (combined)"

Missing files are silently skipped (e.g. if gene-set analysis was not run).
"""

import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ──────────────────────────────────────────────────────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-prefix", required=True,
                    help="Prefix of MAGMA output files (without extension)")
    ap.add_argument("--plots-dir",     required=True,
                    help="Directory to write PNG plots")
    ap.add_argument("--gene-loc",      required=True,
                    help="Path to NCBI*.gene.loc (Entrez ID → symbol map)")
    ap.add_argument("--run-name",      default="",
                    help="Label used in plot titles")
    return ap.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def load_gene_symbols(gene_loc_path: str) -> dict[str, str]:
    """Return {entrez_id: gene_symbol} from NCBI*.gene.loc."""
    sym: dict[str, str] = {}
    try:
        with open(gene_loc_path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 6:
                    sym[parts[0]] = parts[5]
    except FileNotFoundError:
        print(f"WARNING: gene-loc file not found: {gene_loc_path}", file=sys.stderr)
    return sym


def read_genes_out(path: str, gene_sym: dict, max_rows: int | None = None) -> list[dict]:
    """Parse magma *.genes.out file. Returns list of gene dicts sorted by p."""
    rows: list[dict] = []
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 9 or parts[0] in ("GENE", "#"):
                    continue
                try:
                    rows.append({
                        "gene_id": parts[0],
                        "chr":     int(parts[1]),
                        "start":   int(parts[2]),
                        "stop":    int(parts[3]),
                        "nsnps":   int(parts[4]),
                        "zstat":   float(parts[7]),
                        "p":       float(parts[8]),
                        "symbol":  gene_sym.get(parts[0], parts[0]),
                    })
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        return []
    rows.sort(key=lambda r: r["p"])
    return rows[:max_rows] if max_rows else rows


def read_gsa_out(path: str) -> list[dict]:
    """Parse magma *.gsa.out file. Returns list of gene-set dicts sorted by p."""
    rows: list[dict] = []
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#") or line.startswith("VARIABLE"):
                    continue
                parts = line.split()
                if len(parts) < 7:
                    continue
                try:
                    rows.append({
                        "name":   parts[0],
                        "type":   parts[1],
                        "ngenes": int(parts[2]),
                        "beta":   float(parts[3]),
                        "se":     float(parts[5]),
                        "p":      float(parts[6]),
                    })
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        return []
    rows.sort(key=lambda r: r["p"])
    return rows


# ──────────────────────────────────────────────────────────────────────────────
# Plots
# ──────────────────────────────────────────────────────────────────────────────

def plot_gene_manhattan(all_genes: list[dict], run_name: str, out_path: str):
    """Genome-wide gene-based –log10(P) Manhattan plot."""
    if not all_genes:
        print("  Skipping gene Manhattan (no gene data)", file=sys.stderr)
        return

    # Build cumulative chromosome offsets
    chr_sizes: dict[int, int] = {}
    for g in all_genes:
        chr_sizes[g["chr"]] = max(chr_sizes.get(g["chr"], 0), g["stop"])

    chrs = sorted(chr_sizes)
    offset: dict[int, int] = {}
    cum = 0
    for c in chrs:
        offset[c] = cum
        cum += chr_sizes[c] + 5_000_000

    fig, ax = plt.subplots(figsize=(16, 5))
    colors = ["#1f78b4", "#33a02c"]
    for g in all_genes:
        x = offset[g["chr"]] + g["start"]
        y = -math.log10(g["p"]) if g["p"] > 0 else 0
        ax.scatter(x, y, s=6,
                   c=colors[chrs.index(g["chr"]) % 2],
                   alpha=0.5, linewidths=0)

    # Chromosome tick labels
    xticks = [offset[c] + chr_sizes[c] / 2 for c in chrs]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(c) for c in chrs], fontsize=7)
    ax.set_xlabel("Chromosome", fontsize=11)
    ax.set_ylabel("−log₁₀(P)", fontsize=11)
    ax.set_title(f"MAGMA Gene-Based Test — {run_name}", fontsize=13)

    # Thresholds
    n_genes = len(all_genes) or 18_482
    bonf_y = -math.log10(0.05 / n_genes)
    nom_y  = -math.log10(0.05)
    ax.axhline(bonf_y, color="red", linestyle="--", linewidth=0.8,
               label=f"Bonferroni (P={0.05/n_genes:.1e})")
    ax.axhline(nom_y,  color="grey", linestyle=":", linewidth=0.8,
               label="P = 0.05")

    # Annotate top 10 genes
    top10 = sorted(all_genes, key=lambda g: g["p"])[:10]
    for g in top10:
        if g["symbol"]:
            x = offset[g["chr"]] + g["start"]
            y = -math.log10(g["p"]) if g["p"] > 0 else 0
            ax.annotate(g["symbol"], (x, y),
                        textcoords="offset points", xytext=(0, 4),
                        fontsize=7, ha="center", color="black")

    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}", flush=True)


def _clean_name(name: str, prefixes: list[str] = ()) -> str:
    s = name
    for p in prefixes:
        s = s.replace(p, "")
    return s.replace("_", " ").title()


def plot_hallmark_barplot(rows: list[dict], run_name: str, out_path: str):
    if not rows:
        print("  Skipping Hallmark barplot (no data)", file=sys.stderr)
        return
    top = rows[:20]
    names  = [_clean_name(r["name"], ["HALLMARK_"]) for r in top]
    ylogs  = [-math.log10(r["p"]) if r["p"] > 0 else 0 for r in top]
    colors = ["#d73027" if r["beta"] > 0 else "#4575b4" for r in top]
    bonf   = -math.log10(0.05 / max(len(rows), 1))

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(range(len(names)), ylogs[::-1],
            color=colors[::-1], edgecolor="grey", linewidth=0.5)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names[::-1], fontsize=9)
    ax.set_xlabel("−log₁₀(P)", fontsize=11)
    ax.set_title(f"MAGMA Hallmark Gene Sets — {run_name}", fontsize=12)
    ax.axvline(bonf, color="red", linestyle="--", linewidth=1,
               label=f"Bonferroni (P={0.05/max(len(rows),1):.3f})")
    ax.axvline(-math.log10(0.05), color="grey", linestyle=":", linewidth=0.8,
               label="P = 0.05")
    p1 = mpatches.Patch(color="#d73027", label="Positive β (enriched)")
    p2 = mpatches.Patch(color="#4575b4", label="Negative β (depleted)")
    ax.legend(handles=[p1, p2] + ax.get_legend_handles_labels()[0][2:], fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}", flush=True)


def _src_color_c2(name: str) -> str:
    if "KEGG"     in name: return "#e41a1c"
    if "REACTOME" in name: return "#377eb8"
    if name.startswith("WP_") or name.startswith("WP "): return "#4daf4a"
    return "#984ea3"


def plot_c2_bubbleplot(rows: list[dict], run_name: str, out_path: str):
    if not rows:
        print("  Skipping C2 bubbleplot (no data)", file=sys.stderr)
        return
    top  = rows[:25]
    names = [_clean_name(r["name"],
                         ["KEGG_", "REACTOME_", "WP_", "BIOCARTA_", "PID_"])[:55]
             for r in top]
    ylogs  = [-math.log10(r["p"]) if r["p"] > 0 else 0 for r in top]
    sizes  = [max(20, r["ngenes"] * 1.5) for r in top]
    colors = [_src_color_c2(r["name"]) for r in top]
    bonf   = -math.log10(0.05 / max(len(rows), 1))

    fig, ax = plt.subplots(figsize=(7, 10))
    ax.scatter(ylogs, range(len(top)), s=sizes, c=colors,
               alpha=0.8, edgecolors="grey", linewidths=0.5)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("−log₁₀(P)", fontsize=11)
    ax.set_title(f"MAGMA C2 Curated Pathways — {run_name}", fontsize=12)
    ax.axvline(bonf, color="red", linestyle="--", linewidth=0.8,
               label=f"Bonferroni P={0.05/max(len(rows),1):.1e}")
    ax.axvline(-math.log10(0.05), color="grey", linestyle=":", linewidth=0.8,
               label="P = 0.05")
    patches = [
        mpatches.Patch(color="#e41a1c", label="KEGG"),
        mpatches.Patch(color="#377eb8", label="Reactome"),
        mpatches.Patch(color="#4daf4a", label="WikiPathways"),
        mpatches.Patch(color="#984ea3", label="Other"),
    ]
    ax.legend(handles=patches, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}", flush=True)


def plot_c5_barplot(rows: list[dict], run_name: str, out_path: str):
    if not rows:
        print("  Skipping C5 barplot (no data)", file=sys.stderr)
        return
    gobp = [r for r in rows if r["name"].startswith("GOBP_")][:15]
    hp   = [r for r in rows if r["name"].startswith("HP_")][:10]
    top  = sorted(gobp + hp, key=lambda r: r["p"])[:20]
    if not top:
        top = rows[:20]   # fallback: just take top 20 regardless of prefix

    names  = [_clean_name(r["name"], ["GOBP_", "GOCC_", "GOMF_", "HP_"])[:60]
              for r in top]
    ylogs  = [-math.log10(r["p"]) if r["p"] > 0 else 0 for r in top]
    colors = ["#1b9e77" if r["name"].startswith("GOBP_") else "#d95f02"
              for r in top]
    bonf   = -math.log10(0.05 / max(len(rows), 1))

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(range(len(names)), ylogs[::-1],
            color=colors[::-1], edgecolor="grey", linewidth=0.5)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names[::-1], fontsize=9)
    ax.set_xlabel("−log₁₀(P)", fontsize=11)
    ax.set_title(f"MAGMA C5 GO-BP + HPO — {run_name}", fontsize=12)
    ax.axvline(bonf, color="red", linestyle="--", linewidth=0.8,
               label=f"Bonferroni P={0.05/max(len(rows),1):.1e}")
    ax.axvline(-math.log10(0.05), color="grey", linestyle=":", linewidth=0.8,
               label="P = 0.05")
    p1 = mpatches.Patch(color="#1b9e77", label="GO Biological Process")
    p2 = mpatches.Patch(color="#d95f02", label="Human Phenotype Ontology")
    ax.legend(handles=[p1, p2], fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# Summary tables
# ──────────────────────────────────────────────────────────────────────────────

def write_top_genes(genes: list[dict], out_path: str):
    with open(out_path, "w") as fh:
        fh.write("Rank\tEntrezID\tSymbol\tCHR\tNSNPS\tZ\tP\n")
        for i, g in enumerate(genes[:50], 1):
            fh.write(
                f"{i}\t{g['gene_id']}\t{g['symbol']}\t{g['chr']}\t"
                f"{g['nsnps']}\t{g['zstat']:.4f}\t{g['p']:.4e}\n"
            )
    print(f"  Saved: {out_path}", flush=True)


def write_top_genesets(rows: list[dict], label: str, out_path: str):
    with open(out_path, "w") as fh:
        fh.write("Rank\tGeneSet\tNGenes\tBeta\tSE\tP\n")
        for i, r in enumerate(rows[:50], 1):
            fh.write(
                f"{i}\t{r['name']}\t{r['ngenes']}\t"
                f"{r['beta']:.4f}\t{r['se']:.4f}\t{r['p']:.4e}\n"
            )
    print(f"  Saved: {out_path}", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.plots_dir, exist_ok=True)

    print(f"Loading gene symbols from {args.gene_loc} …", flush=True)
    gene_sym = load_gene_symbols(args.gene_loc)
    print(f"  {len(gene_sym):,} gene symbols loaded", flush=True)

    prefix = args.output_prefix
    run_name = args.run_name

    # ── 1. Gene Manhattan ─────────────────────────────────────────────────────
    genes_out_path = f"{prefix}.genes.out"
    print(f"\n[1/5] Gene-based Manhattan", flush=True)
    all_genes = read_genes_out(genes_out_path, gene_sym)
    if all_genes:
        print(f"  {len(all_genes):,} genes loaded", flush=True)
        plot_gene_manhattan(
            all_genes, run_name,
            os.path.join(args.plots_dir, "gene_manhattan.png"),
        )
        write_top_genes(
            all_genes,
            os.path.join(os.path.dirname(prefix), "top_genes.tsv"),
        )
    else:
        print(f"  {genes_out_path} not found or empty — skipping", flush=True)

    # ── 2-4. Gene-set plots ───────────────────────────────────────────────────
    gsa_configs = [
        ("h",  "Hallmark",     "hallmark_barplot.png",    "top_genesets_Hallmark.tsv"),
        ("c2", "C2_Curated",   "c2_bubbleplot.png",       "top_genesets_C2_Curated.tsv"),
        ("c5", "C5_GO_HPO",    "c5_barplot.png",          "top_genesets_C5_GO_HPO.tsv"),
    ]
    for step_i, (coll, label, plot_fname, tsv_fname) in enumerate(gsa_configs, 2):
        gsa_path = f"{prefix}_{coll}.gsa.out"
        print(f"\n[{step_i}/5] {label} gene-set analysis", flush=True)
        rows = read_gsa_out(gsa_path)
        if not rows:
            print(f"  {gsa_path} not found or empty — skipping", flush=True)
            continue
        print(f"  {len(rows):,} gene sets loaded", flush=True)

        plot_path = os.path.join(args.plots_dir, plot_fname)
        tsv_path  = os.path.join(os.path.dirname(prefix), tsv_fname)

        if coll == "h":
            plot_hallmark_barplot(rows, run_name, plot_path)
        elif coll == "c2":
            plot_c2_bubbleplot(rows, run_name, plot_path)
        else:
            plot_c5_barplot(rows, run_name, plot_path)

        write_top_genesets(rows, label, tsv_path)

    print("\nAll done.", flush=True)


if __name__ == "__main__":
    main()
