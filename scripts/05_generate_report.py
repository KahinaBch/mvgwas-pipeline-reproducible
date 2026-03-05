#!/usr/bin/env python3
# =============================================================================
# scripts/05_generate_report.py  —  Step V: Comprehensive run report
# =============================================================================
"""
Reads pipeline outputs, run_stats.tsv, and log files to generate a detailed
Markdown (+ optional HTML) report summarising every step of the mvgwas run.
"""
import argparse
import csv
import datetime
import os
import re
import sys
from pathlib import Path

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate mvgwas-pipeline run report")
parser.add_argument("--output-dir",   required=True)
parser.add_argument("--run-name",     required=True)
parser.add_argument("--genome-build", default="GRCh38")
parser.add_argument("--pval-col",     default="P_manta")
parser.add_argument("--gw-thresh",    default="5e-8")
parser.add_argument("--sug-thresh",   default="1e-5")
parser.add_argument("--config-file",  default=None)
parser.add_argument("--log-file",     default=None)
parser.add_argument("--top-n",        default=20, type=int,
                    help="Number of top SNPs to include in report table")
args = parser.parse_args()

OUT      = Path(args.output_dir)
RUN      = args.run_name
BUILD    = args.genome_build
PVCOL    = args.pval_col
GW       = args.gw_thresh
SUG      = args.sug_thresh
TOP_N    = args.top_n
NOW      = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

REPORT_MD   = OUT / f"{RUN}_report.md"
REPORT_HTML = OUT / f"{RUN}_report.html"

# =============================================================================
# Helper functions
# =============================================================================

def read_metadata(meta_path):
    """Read KEY=VALUE metadata file → dict."""
    d = {}
    if not meta_path.exists():
        return d
    for line in meta_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            k, _, v = line.partition("=")
            d[k.strip()] = v.strip()
    return d


def read_stats(stats_path):
    """Read run_stats.tsv → dict (last value wins per key)."""
    d = {}
    if not stats_path.exists():
        return d
    with open(stats_path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2:
                d[row[0]] = row[1]
    return d


def read_report_file(path):
    """Read a text report file. Return placeholder if missing."""
    if path and Path(path).exists():
        return Path(path).read_text().strip()
    return "_File not found._"


def tail_log(path, n=100):
    """Return last n lines of a log file."""
    if not path or not Path(path).exists():
        return "_Log file not found._"
    lines = Path(path).read_text().splitlines()
    tail = lines[-n:]
    return "\n".join(tail)


def read_top_snps(tsv_path, n):
    """Return list of dicts (first n rows of TSV)."""
    if not tsv_path or not Path(tsv_path).exists():
        return [], []
    rows = []
    header = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        for i, row in enumerate(reader):
            if i >= n:
                break
            rows.append(row)
    return header, rows


def md_table(header, rows):
    """Render list-of-dicts as Markdown table (first 8 cols max)."""
    if not header or not rows:
        return "_No data available._"
    cols = header[:8]
    lines = []
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("|" + "|".join([" --- " for _ in cols]) + "|")
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(c, "")) for c in cols) + " |")
    return "\n".join(lines)


def img_ref(path, alt="plot", width=700):
    """Render an image reference if file exists."""
    if path and Path(path).exists():
        return f'<img src="{path}" alt="{alt}" width="{width}"/>\n'
    return f'_Image not found: {path}_\n'


def section_status(condition, label):
    return f"✅ {label}" if condition else f"⚠️  {label} (not found)"


# =============================================================================
# Load all data sources
# =============================================================================
meta        = read_metadata(OUT / f"{RUN}_run_metadata.txt")
stats       = read_stats(OUT / "run_stats.tsv")
qc_dir      = OUT / "qc"
plots_dir   = meta.get("PLOTS_DIR", str(OUT / "plots"))

val_report  = read_report_file(meta.get("VALIDATION_REPORT",
                str(OUT / "inputs" / "input_validation_report.txt")))
merge_qc    = read_report_file(meta.get("MERGE_QC_REPORT",
                str(qc_dir / "merge_qc_report.txt")))

merged_results  = meta.get("MERGED_RESULTS", "")
annotated_snps  = meta.get("ANNOTATED_SNPS", "")
top_snps_path   = annotated_snps or meta.get("TOP_SNPS", "")
plot_mht        = meta.get("PLOT_MANHATTAN", str(Path(plots_dir) / f"manhattan_{RUN}.png"))
plot_qq         = meta.get("PLOT_QQ",        str(Path(plots_dir) / f"qq_{RUN}.png"))
plot_reg        = meta.get("PLOT_REGIONAL",  str(Path(plots_dir) / f"regional_{RUN}.png"))

header_top, rows_top = read_top_snps(top_snps_path, TOP_N)

# Key statistics
n_samples     = meta.get("SAMPLES_INTERSECTION", stats.get("samples_intersection", "?"))
n_total_assoc = meta.get("TOTAL_ASSOCIATIONS",   stats.get("total_associations", "?"))
n_top_extract = stats.get("top_n_extracted", "?")
chrs_present  = stats.get("chromosomes_present", "?")
chrs_missing  = stats.get("chromosomes_missing",  "0")

# Config file excerpt
config_text = ""
if args.config_file and Path(args.config_file).exists():
    config_text = Path(args.config_file).read_text().strip()
elif (OUT / f"{RUN}.conf").exists():
    config_text = (OUT / f"{RUN}.conf").read_text().strip()

# Determine log path
log_path = args.log_file or meta.get("PIPELINE_LOG",
           str(OUT / "logs" / f"{RUN}_pipeline.log"))
log_tail = tail_log(log_path, 150)

# =============================================================================
# Build Markdown
# =============================================================================
md_parts = []

def H(n, text): return "#" * n + " " + text
def HR(): return "\n---\n"


# ── Title & Header ─────────────────────────────────────────────────────────────
md_parts.append(f"""# mvgwas-pipeline Run Report

**Run name:** `{RUN}`  
**Genome build:** {BUILD}  
**Report generated:** {NOW}  
**Pipeline output:** `{OUT}`
""")

md_parts.append(HR())

# ── Executive Summary ─────────────────────────────────────────────────────────
md_parts.append(H(2, "Executive Summary"))

gw_count = "?"
if rows_top and PVCOL in header_top:
    try:
        gw_count = sum(1 for r in rows_top if float(r.get(PVCOL, 1)) < float(GW))
    except Exception:
        pass

md_parts.append(f"""
| Parameter | Value |
| --- | --- |
| Sample size (intersection) | {n_samples} |
| Total associations tested | {n_total_assoc} |
| Chromosomes with results | {chrs_present} / 22 |
| Missing chromosomes | {chrs_missing} |
| Top SNPs extracted | {n_top_extract} |
| GW-significant (p < {GW}) | {gw_count} |
| P-value column | `{PVCOL}` |
| Significance threshold | {GW} |
| Suggestive threshold | {SUG} |
""")

md_parts.append(HR())

# ── Step I: Environment ───────────────────────────────────────────────────────
md_parts.append(H(2, "Step I — Environment & Setup"))
step1_log = tail_log(str(OUT / "logs" / "01_setup.log"), 60)
md_parts.append(f"""
```
{step1_log}
```
""")
md_parts.append(HR())

# ── Step II: Input QC ─────────────────────────────────────────────────────────
md_parts.append(H(2, "Step II — Input Validation & Sample QC"))
md_parts.append(f"""
{val_report}
""")
md_parts.append(HR())

# ── Step III: Pipeline execution ─────────────────────────────────────────────
md_parts.append(H(2, "Step III — Pipeline Execution"))
md_parts.append(f"""
**Merged results file:** `{merged_results}`  
**Total associations:** {n_total_assoc}  
**Chromosomes present:** {chrs_present}  
**Chromosomes missing:** {chrs_missing}  

### Merge QC Report

```
{merge_qc}
```
""")
md_parts.append(HR())

# ── Step IV: Post-processing ──────────────────────────────────────────────────
md_parts.append(H(2, "Step IV — Top SNPs, rsID Annotation & Plots"))

md_parts.append(f"""
**Top SNPs file:** `{top_snps_path}`  
**rsID annotation:** {'✅ done' if annotated_snps else '⚠️ not performed'}  

### Top {min(TOP_N, len(rows_top))} Associations

{md_table(header_top, rows_top)}
""")

md_parts.append(H(3, "Manhattan Plot"))
md_parts.append(img_ref(plot_mht, "Manhattan", 900))

md_parts.append(H(3, "QQ Plot"))
md_parts.append(img_ref(plot_qq, "QQ plot", 600))

md_parts.append(H(3, "Regional Manhattan (Top Locus)"))
md_parts.append(img_ref(plot_reg, "Regional Manhattan", 800))

md_parts.append(HR())

# ── Step V: Report metadata ───────────────────────────────────────────────────
md_parts.append(H(2, "Step V — Report & Statistics"))

if stats:
    stat_lines = ["| Key | Value |", "| --- | --- |"]
    for k, v in stats.items():
        stat_lines.append(f"| {k} | {v} |")
    md_parts.append("\n".join(stat_lines) + "\n")
else:
    md_parts.append("_run_stats.tsv not found._\n")

md_parts.append(HR())

# ── Config used ───────────────────────────────────────────────────────────────
if config_text:
    md_parts.append(H(2, "Run Configuration"))
    md_parts.append(f"```bash\n{config_text}\n```\n")
    md_parts.append(HR())

# ── Run metadata dump ─────────────────────────────────────────────────────────
if meta:
    md_parts.append(H(2, "Run Metadata"))
    meta_lines = ["| Key | Value |", "| --- | --- |"]
    for k, v in meta.items():
        meta_lines.append(f"| `{k}` | `{v}` |")
    md_parts.append("\n".join(meta_lines) + "\n")
    md_parts.append(HR())

# ── Pipeline log tail ─────────────────────────────────────────────────────────
md_parts.append(H(2, "Pipeline Log (last 150 lines)"))
md_parts.append(f"```\n{log_tail}\n```\n")

# =============================================================================
# Write Markdown
# =============================================================================
full_md = "\n".join(md_parts)
REPORT_MD.write_text(full_md)
print(f"[OK] Markdown report: {REPORT_MD}")

# =============================================================================
# Convert to HTML (inline CSS, no external dependencies)
# =============================================================================
try:
    import subprocess
    # Try pandoc first
    result = subprocess.run(
        ["pandoc", "-f", "markdown", "-t", "html",
         "--standalone", "--embed-resources",
         "-o", str(REPORT_HTML), str(REPORT_MD)],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        print(f"[OK] HTML report (pandoc): {REPORT_HTML}")
    else:
        raise RuntimeError("pandoc failed")
except Exception:
    # Minimal fallback using Python stdlib
    html_head = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>mvgwas Report — {RUN}</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1100px;margin:40px auto;padding:0 20px;color:#222}}
h1{{color:#1a1a6e}} h2{{color:#2c2c8e;border-bottom:2px solid #dde}} h3{{color:#444}}
table{{border-collapse:collapse;width:100%;margin:12px 0}}
th,td{{border:1px solid #ccc;padding:6px 10px;font-size:12px}}
th{{background:#f0f0f0}} tr:nth-child(even){{background:#fafafa}}
code,pre{{background:#f5f5f5;padding:2px 6px;border-radius:3px;font-size:12px}}
pre{{padding:12px;overflow-x:auto;max-height:500px;border:1px solid #ddd}}
img{{max-width:100%;border:1px solid #ddd;margin:8px 0}}
hr{{border:none;border-top:2px solid #eee;margin:30px 0}}
</style></head><body>
"""
    # Simple markdown → HTML (headers, code blocks, tables, horizontal rules)
    html_body = full_md
    html_body = re.sub(r'^### (.+)$', r'<h3>\1</h3>', html_body, flags=re.M)
    html_body = re.sub(r'^## (.+)$',  r'<h2>\1</h2>', html_body, flags=re.M)
    html_body = re.sub(r'^# (.+)$',   r'<h1>\1</h1>', html_body, flags=re.M)
    html_body = re.sub(r'^---$', r'<hr>', html_body, flags=re.M)
    html_body = re.sub(r'```[a-z]*\n(.*?)```', lambda m:
        f'<pre><code>{m.group(1)}</code></pre>', html_body, flags=re.S)
    html_body = re.sub(r'`([^`]+)`', r'<code>\1</code>', html_body)
    html_body = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', html_body)
    html_body = html_body.replace("\n\n", "</p><p>")
    REPORT_HTML.write_text(html_head + "<p>" + html_body + "</p></body></html>")
    print(f"[OK] HTML report (fallback): {REPORT_HTML}")

print(f"[OK] Step V complete. Reports in: {OUT}")
