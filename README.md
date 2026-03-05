# mvgwas-pipeline

An automated, reproducible pipeline wrapper for
[KahinaBch/mvgwas-nf](https://github.com/KahinaBch/mvgwas-nf) (branch `kb/dsl2-conversion`).

Performs five sequential steps: environment setup → input QC → chromosome-parallel GWAS
execution → post-processing (top SNPs, rsID, Manhattan/QQ/Regional plots) → report generation.

---
## Quick overview

### Inputs

Three files are required, all specified in your config:

| File | Config key | Format | Notes |
|------|-----------|--------|-------|
| Genotype VCF | `VCF_FILE` | bgzipped VCF (`.vcf.gz`) + tabix index (`.vcf.gz.tbi`) | Multi-sample; one variant per row; sample IDs must match phenotype/covariate files |
| Phenotype file | `PHENOTYPE_FILE` | Tab-separated (`.tsv` or `.tsv.gz`), **header required** | First column = sample ID; remaining columns = phenotype traits (one per column); no row index |
| Covariate file | `COVARIATE_FILE` | Tab-separated (`.tsv` or `.tsv.gz`), **header required** | First column = sample ID; remaining columns = covariates (age, sex, PCs, …); same sample set as phenotype file |

**Sample ID matching:** the pipeline computes the three-way intersection of sample IDs found in the VCF header, the phenotype file, and the covariate file. Only samples present in all three are retained for the analysis. The pipeline aborts if fewer than 50 samples pass this filter.

**Chromosome naming:** both `chr1` and `1` styles are auto-detected from the VCF header.

---

### Outputs

All outputs are written under `OUTPUT_DIR` (set in config):

```
OUTPUT_DIR/
├── inputs/                          ← Step II
│   ├── samples_vcf.txt              # All sample IDs extracted from the VCF
│   ├── samples_intersection.txt     # Final sample list (three-way overlap)
│   ├── phenotype_filtered.tsv       # Phenotype file restricted to intersection samples
│   ├── covariate_filtered.tsv       # Covariate file restricted to intersection samples
│   └── input_validation_report.txt  # Detailed QC log (counts, duplicates, NAs)
│
├── chr_vcfs/                        ← Step III (intermediate)
│   ├── chr1.vcf.gz + .tbi           # Per-chromosome VCF subsets
│   └── ...
│
├── chr_results/                     ← Step III (intermediate)
│   ├── chr1/result/mvgwas_chr1.tsv  # Raw mvgwas-nf output per chromosome
│   └── ...
│
├── results/                         ← Steps III & IV
│   ├── mvgwas_merged.tsv            # Full merged association results (all chromosomes)
│   ├── top_1000_snps.tsv            # Top N associations sorted by p-value
│   └── top_1000_snps_rsid.tsv       # Same, with rsID annotation column added
│
├── qc/                              ← Step III
│   └── merge_qc_report.txt          # Per-chromosome variant counts, missing chr list
│
├── plots/                           ← Step IV
│   ├── manhattan_<RUN_NAME>.png     # Genome-wide Manhattan plot
│   ├── qq_<RUN_NAME>.png            # QQ plot with λ GC
│   └── regional_<RUN_NAME>.png      # Regional Manhattan ± window around top locus
│
├── logs/                            ← All steps
│   ├── <RUN_NAME>_pipeline.log      # Master timestamped log
│   ├── chr1.log … chr22.log         # Per-chromosome nextflow logs
│   └── plotting.log                 # R plotting log
│
├── <RUN_NAME>_run_metadata.txt      # Key=value store of all runtime paths & stats
├── run_stats.tsv                    # Tab-separated summary statistics table
├── <RUN_NAME>_report.md             ← Step V — human-readable run report
└── <RUN_NAME>_report.html           ← Step V — same report as standalone HTML
```

---

## Requirements

| Tool | Min version | Notes |
|------|-------------|-------|
| Bash | 4.0+ | macOS ships with 3.x — install via Homebrew |
| Java | 11 (or 8) | Java 8 triggers auto NXF_VER=22.04.0 |
| Nextflow | 21.04+ | auto-installed on first run |
| bcftools | 1.12+ | |
| bgzip / tabix | same | shipped with htslib |
| Rscript | 4.0+ | R packages auto-installed |
| python3 | 3.8+ | pip packages auto-installed |
| Singularity | 3.x | optional; or Docker; or none |

---

## Quick Start

```bash
# 1. Clone this repository
git clone https://github.com/YOUR_ORG/mvgwas-pipeline.git
cd mvgwas-pipeline

# 2. Copy and fill in the config template
cp config.conf.template my_run.conf
# Edit my_run.conf — set VCF_FILE, PHENOTYPE_FILE, COVARIATE_FILE, OUTPUT_DIR, PIPELINE_DIR

# 3. Launch
bash run.sh --config my_run.conf
```

---

## Command-Line Options

```
Usage: bash run.sh [OPTIONS]

Options:
  --config <file>          Path to configuration file (required)
  --skip-test              Skip Step I (environment test + pipeline self-test)
  --build <GRCh37|GRCh38>  Override genome build (also set in config)
  --steps <1,2,3,4,5>      Run only specific steps (e.g. --steps 4,5)
  --resume                 Skip steps whose checkpoint already exists
  --dry-run                Print commands without executing them
  --verbose                Enable DEBUG-level logging
  --help                   Print this help message
```

---

## Five Pipeline Steps

### Step I — Environment & Self-Test
- Checks all required tools (java, nextflow, bcftools, bgzip, tabix, Rscript, python3)
- Auto-detects Java 8 and pins Nextflow to 22.04.0 for compatibility
- Verifies that `PIPELINE_DIR` is on the correct git branch (`kb/dsl2-conversion`)
- Checks and auto-installs required R packages (data.table, ggplot2, ggrepel, scales, cowplot, jsonlite, optparse)
- Checks and auto-installs required Python packages (pandas, numpy, scipy, matplotlib, jinja2)
- Runs a pipeline self-test using the bundled test data in `mvgwas-nf/data/`
- **Skip with `--skip-test`** or by setting `SKIP_TEST=true` in your config

### Step II — Input Validation & Sample Overlap
- Validates VCF (with index), phenotype TSV, and covariate TSV
- Extracts VCF sample IDs using `bcftools query -l`
- Detects chromosome prefix (`chr1` vs `1`)
- Checks for duplicate IDs, missing values, and column completeness
- Computes the three-way intersection of samples across VCF × phenotype × covariates
- Filters phenotype and covariate files to the intersection set
- Reports exclusion counts per source
- Aborts if the sample intersection is < 50

### Step III — Chromosome-Parallel GWAS Execution
- Splits the VCF by chromosome using `bcftools view -r`
- Runs `nextflow run mvgwas.nf` per chromosome with the following parameters:
  `--geno`, `--pheno`, `--cov`, `--dir`, `--out`, `--l` (chunk size), `--t` (transformation), `--i` (interaction), `--ng` (min individuals)
- **Local mode:** runs up to `LOCAL_PARALLEL` (default: 4) chromosomes in parallel
- **SLURM mode:** submits a job array (1–22) and a dependent merge job
- Merges all chromosome results (header once + data rows)
- Performs merge QC: chromosome completeness, variant counts, column validation

### Step IV — Post-Processing
- Extracts the top `TOP_N_SNPS` (default: 1000) associations sorted by `PVAL_COL`
- Annotates rsIDs using `bcftools annotate` against a local or remote dbSNP VCF
  - **GRCh38:** `GCF_000001405.40.gz` (NCBI FTP)
  - **GRCh37:** `GCF_000001405.25.gz` (NCBI FTP)
  - Set `DBSNP_VCF` in your config to use a local file (much faster)
- Generates:
  - **Manhattan plot** — full genome, coloured by chromosome, labelled top SNPs
  - **QQ plot** — observed vs. expected −log₁₀(P) with lambda GC
  - **Regional Manhattan** — ±`REGIONAL_WINDOW` bp around the top locus, with UCSC gene track

### Step V — Report Generation
- Aggregates all intermediate results, QC reports, and statistics
- Produces a comprehensive **Markdown** report (+ **HTML** via pandoc or built-in fallback)
- Sections: executive summary, run config, step-by-step QC details, top loci table, embedded plots, full pipeline log tail

---

## Configuration Reference

Copy `config.conf.template` and edit:

```bash
# === Run identity ===
RUN_NAME="my_analysis"
OUTPUT_DIR="/path/to/output"

# === Inputs ===
VCF_FILE="/path/to/genotypes.vcf.gz"
PHENOTYPE_FILE="/path/to/phenotypes.tsv"
COVARIATE_FILE="/path/to/covariates.tsv"

# === Pipeline location ===
PIPELINE_DIR="/path/to/mvgwas-nf"
PIPELINE_BRANCH="kb/dsl2-conversion"  # default

# === Analysis parameters ===
CHUNK_SIZE=500
TRANSFORMATION=none   # none | rank | log | sqrt
INTERACTION=none
MIN_INDIVIDUALS=10
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

# === Container ===
CONTAINER=singularity  # singularity | docker | none

# === Genome & annotation ===
GENOME_BUILD=GRCh38   # GRCh37 | GRCh38
DBSNP_VCF=""          # optional: path to local dbSNP VCF.gz

# === Post-processing ===
TOP_N_SNPS=1000
TOP_LABEL_N=20
PVAL_COL=P_manta
GW_THRESH=5e-8
SUG_THRESH=1e-5
REGIONAL_WINDOW=2000000

# === Execution ===
USE_SLURM=false
LOCAL_PARALLEL=4
```

For SLURM, configure the per-chromosome job and the merge job separately:

```bash
# === SLURM — per-chromosome array job ===
SLURM_PARTITION="genoa"       # partition / queue name on your cluster
SLURM_TIME="24:00:00"         # wall-time per chromosome job
SLURM_MEM="32G"               # memory per chromosome job
SLURM_CPUS=4                  # CPUs per chromosome job
SLURM_MAIL=""                 # email for failure notifications (empty = off)
SLURM_EXTRA_ARGS=""           # any extra #SBATCH directives (e.g. --account=...)

# === SLURM — merge job (submitted with afterok dependency) ===
SLURM_MERGE_TIME="02:00:00"   # wall-time for the merge step
SLURM_MERGE_MEM="16G"         # memory for the merge step
SLURM_MERGE_CPUS=2            # CPUs for the merge step
```

The chromosome job array is built automatically from the `CHROMOSOMES` variable, so setting e.g. `CHROMOSOMES="1 2 3"` will submit `--array=1,2,3` rather than the full 1–22.

---

## Output Structure

```
OUTPUT_DIR/
├── inputs/
│   ├── samples_vcf.txt
│   ├── samples_intersection.txt
│   ├── phenotype_filtered.tsv
│   ├── covariate_filtered.tsv
│   └── input_validation_report.txt
├── chr_vcfs/
│   ├── chr1.vcf.gz  (.tbi)
│   └── ...
├── chr_results/
│   ├── chr1/result/mvgwas_chr1.tsv
│   └── ...
├── results/
│   ├── mvgwas_merged.tsv
│   ├── top_1000_snps.tsv
│   └── top_1000_snps_rsid.tsv
├── qc/
│   └── merge_qc_report.txt
├── plots/
│   ├── manhattan_<RUN_NAME>.png
│   ├── qq_<RUN_NAME>.png
│   └── regional_<RUN_NAME>.png
├── logs/
│   ├── <RUN_NAME>_pipeline.log
│   ├── chr1.log … chr22.log
│   └── plotting.log
├── <RUN_NAME>_run_metadata.txt
├── run_stats.tsv
├── <RUN_NAME>_report.md
└── <RUN_NAME>_report.html
```

---

## Resuming a Run

If a step was interrupted:

```bash
bash run.sh --config my_run.conf --resume
```

`--resume` skips any step whose checkpoint file already exists under `OUTPUT_DIR/.checkpoints/`.

To re-run only post-processing and the report:

```bash
bash run.sh --config my_run.conf --steps 4,5 --resume
```

---

## Input File Formats

### Phenotype / Covariate TSV
- Tab-separated, **header row required**
- First column must be **sample ID** (matching VCF sample names)
- All phenotype and covariate columns are passed directly to `mvgwas-nf`

### VCF
- Bgzipped + tabix-indexed (`.vcf.gz` + `.vcf.gz.tbi`)
- Sample IDs must overlap with phenotype / covariate files

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Nextflow not found` | Install: `curl -s https://get.nextflow.io | bash` then move to PATH |
| `bcftools: command not found` | `conda install -c bioconda bcftools` or `brew install bcftools` |
| Java 8 detected | Auto-handled (NXF_VER=22.04.0); or upgrade Java |
| SLURM jobs queued forever | Check `SLURM_PARTITION` in config |
| rsID annotation slow | Set `DBSNP_VCF` to a local pre-downloaded dbSNP file |
| Missing chromosome result | Check `logs/chrN.log` for that chromosome |
| Regional plot missing gene track | UCSC API may be unavailable; plot is still generated without genes |

---

## Citation

If you use this pipeline, please cite:

> [mvgwas-nf: multi-variate GWAS Nextflow pipeline](https://github.com/dgarrimar/mvgwas-nf)  
> Original pipeline by Diego Garrido-Martín.  
> Fork branch `kb/dsl2-conversion` by KahinaBch.

---

## License

MIT License — see [LICENSE](LICENSE).
