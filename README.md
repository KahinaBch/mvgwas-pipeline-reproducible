# mvgwas-pipeline

An automated, reproducible pipeline wrapper for
[KahinaBch/mvgwas-nf](https://github.com/KahinaBch/mvgwas-nf) (branch `kb/dsl2-conversion`).

Performs six sequential steps: environment setup → input QC → chromosome-parallel GWAS
execution → post-processing (top SNPs, rsID, Manhattan/QQ/Regional plots) →
enrichment analysis (GWAS Catalog mapping + MAGMA gene/gene-set test) → report generation.

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
├── chr_results/                         ← Step III — combined stratum
│   ├── chr1/result/mvgwas_chr1.tsv      # Raw mvgwas-nf output per chromosome
│   └── ...
├── chr_results_male/                    ← Step III — male stratum (if SEX_STRATIFIED)
│   └── chr1/result/mvgwas_chr1.tsv
├── chr_results_female/                  ← Step III — female stratum (if SEX_STRATIFIED)
│   └── ...
│
├── results/                             ← Steps III & IV — combined stratum
│   ├── mvgwas_merged.tsv                # Full merged association results
│   ├── top_1000_snps.tsv                # Top N associations sorted by p-value
│   └── top_1000_snps_rsid.tsv           # Same, with rsID annotation column
├── results_male/                        ← Steps III & IV — male stratum (if SEX_STRATIFIED)
│   ├── mvgwas_merged.tsv
│   ├── top_1000_snps.tsv
│   └── top_1000_snps_rsid.tsv
├── results_female/                      ← Steps III & IV — female stratum (if SEX_STRATIFIED)
│   └── ...
│
├── qc/                                  ← Step III
│   └── merge_qc_report.txt              # Per-chromosome variant counts, missing chr list
│
├── plots/                               ← Step IV — combined stratum
│   ├── manhattan_<RUN_NAME>.png         # Genome-wide Manhattan plot
│   ├── qq_<RUN_NAME>.png                # QQ plot with λ GC
│   └── regional_<RUN_NAME>.png          # Regional Manhattan ± window around top locus
├── plots_male/                          ← Step IV — male stratum (if SEX_STRATIFIED)
│   └── manhattan_<RUN_NAME>_male.png  …
├── plots_female/                        ← Step IV — female stratum (if SEX_STRATIFIED)
│   └── manhattan_<RUN_NAME>_female.png …
│
├── logs/                            ← All steps
│   ├── <RUN_NAME>_pipeline.log      # Master timestamped log
│   ├── chr1.log … chr22.log         # Per-chromosome nextflow logs
│   └── plotting.log                 # R plotting log
│
├── <RUN_NAME>_run_metadata.txt      # Key=value store of all runtime paths & stats
├── run_stats.tsv                    # Tab-separated summary statistics table
├── gwas_catalog/                    ← Step V-A — GWAS Catalog mapping
│   ├── combined/
│   │   ├── gwas_catalog_combined.md   # Markdown report per locus
│   │   └── gwas_catalog_combined.tsv  # TSV summary (lead SNPs + catalog hits)
│   ├── male/   (if SEX_STRATIFIED)
│   └── female/ (if SEX_STRATIFIED)
│
├── magma/                           ← Step V-B — MAGMA enrichment
│   ├── combined/
│   │   ├── input/snp_loc.txt         # SNP CHR BP for --annotate
│   │   ├── input/pval.txt            # SNP P for --pval
│   │   ├── output/<RUN>_combined.genes.out   # gene-based results
│   │   ├── output/<RUN>_combined_h.gsa.out   # Hallmark gene-set results
│   │   ├── output/<RUN>_combined_c2.gsa.out  # C2 curated pathway results
│   │   ├── output/<RUN>_combined_c5.gsa.out  # C5 GO+HPO results
│   │   ├── plots/gene_manhattan.png
│   │   ├── plots/hallmark_barplot.png
│   │   ├── plots/c2_bubbleplot.png
│   │   ├── plots/c5_barplot.png
│   │   └── top_genes.tsv / top_genesets_*.tsv
│   ├── male/   (if SEX_STRATIFIED)
│   └── female/ (if SEX_STRATIFIED)
│
├── <RUN_NAME>_report.md             ← Step VI — human-readable run report
└── <RUN_NAME>_report.html           ← Step VI — same report as standalone HTML
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
  --sex-stratified         Run sex-stratified analysis (male + female + combined)
  --help                   Print this help message
```

---

## Six Pipeline Steps

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
- **Sex stratification** (when `SEX_STRATIFIED=true` or `--sex-stratified`):
  - Auto-detects the sex/gender column in the covariate file (looks for `sex`, `gender`, case-insensitive)
  - Auto-detects coding: `1/2`, `M/F`, `male/female`, or `0/1`; overrideable via `SEX_MALE_CODE` / `SEX_FEMALE_CODE`
  - Produces separate filtered phenotype + covariate files for males and females
  - Writes `samples_male.txt` and `samples_female.txt` to `inputs/`
  - Warns if either stratum has < 50 samples

### Step III — Chromosome-Parallel GWAS Execution
- Splits the VCF by chromosome using `bcftools view -r`
- Runs `nextflow run mvgwas.nf` per chromosome with the following parameters:
  `--geno`, `--pheno`, `--cov`, `--dir`, `--out`, `--l` (chunk size), `--t` (transformation), `--i` (interaction), `--ng` (min individuals)
- **Local mode:** runs up to `LOCAL_PARALLEL` (default: 4) chromosomes in parallel
- **SLURM mode:** submits, per stratum, a job array (1–22) and a dependent merge job, plus one job that
  runs steps IV–VI after all merges succeed. `run.sh` exits after submitting; no second call is needed
  (the submitted job IDs are in `OUTPUT_DIR/.slurm_pending` until step IV starts)
- Each chromosome job starts Nextflow in its own folder (`chr_results*/chrN/`), so concurrent runs do not clash
- mvgwas-nf tests a variant only if all three genotype groups have ≥ `MIN_INDIVIDUALS` samples
  (and ≤ 5 % missing genotypes) — in practice a MAF filter of roughly √(`MIN_INDIVIDUALS`/N)
- Merges all chromosome results (header once + data rows)
- Performs merge QC: chromosome completeness, variant counts, column validation

### Step IV — Post-Processing
- Extracts the top `TOP_N_SNPS` (default: 1000) associations sorted by `PVAL_COL`
- Annotates rsIDs using `bcftools annotate` against a local or remote dbSNP VCF
  - **GRCh38:** `GCF_000001405.40.gz` (NCBI FTP)
  - **GRCh37:** `GCF_000001405.25.gz` (NCBI FTP)
  - Set `DBSNP_VCF` in your config to use a local file (much faster, and required where
    the node running step IV has no internet). Contigs may be `1..22` or NCBI-style (`NC_000001.11`)
  - `scripts/make_dbsnp_subset.sh` builds a small local file with only the positions step IV needs
- Generates:
  - **Manhattan plot** — full genome, coloured by chromosome, labelled top SNPs
  - **QQ plot** — observed vs. expected −log₁₀(P) with lambda GC
  - **Regional Manhattan** — ±`REGIONAL_WINDOW` bp around the top locus, with UCSC gene track

### Step V — Enrichment Analysis *(on by default; skip with `--skip-enrichment`)*

#### 5A — GWAS Catalog mapping
- Clusters top SNPs into loci (configurable clumping window, default 500 kb per chr)
- Picks the lead SNP (lowest p-value) per locus and queries the [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas) REST API:
  - **Direct SNP associations** — whether the exact rsID has ever been reported
  - **Region context** — any catalog variants within ±500 kb (biological context)
- Writes a per-stratum Markdown report + TSV summary to `gwas_catalog/<stratum>/`
- Handles offline / rate-limited scenarios gracefully (writes a placeholder report)

#### 5B — MAGMA gene & gene-set enrichment
- Prepares MAGMA input files (SNP-location + p-value) from the merged results
- Annotates SNPs → genes (`magma --annotate`, default ±10 kb window)
- Runs the gene-based association test (`snp-wise=mean` model)
- Optionally runs gene-set analyses against three MSigDB collections
  (when `MAGMA_GENESETS_DIR` is set):
  - **Hallmark** (50 gene sets)
  - **C2** curated pathways (KEGG, Reactome, WikiPathways …)
  - **C5** GO Biological Process + Human Phenotype Ontology
- Produces four plots and summary TSVs in `magma/<stratum>/plots/`

### Step VI — Report Generation
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

# === Sex-stratified analysis ===
SEX_STRATIFIED=false          # true to run male + female + combined
SEX_COL=""                    # sex column in covariate file (auto-detected if empty)
SEX_MALE_CODE=""              # e.g. 1, M, male  (auto-detected if empty)
SEX_FEMALE_CODE=""            # e.g. 2, F, female (auto-detected if empty)

# === Enrichment analysis (Step V) ===
SKIP_ENRICHMENT=false         # true to skip GWAS Catalog + MAGMA entirely
GWAS_CATALOG_TOP_N=20         # top lead SNPs to query
MAGMA_PATH=""                 # path to magma binary (auto-detect if empty)
MAGMA_REF_DIR=""              # dir with NCBI38.gene.loc + g1000_eur.*
MAGMA_GENESETS_DIR=""         # dir with MSigDB .gmt files (skip gene-sets if empty)
MAGMA_SAMPLE_N=""             # sample N for MAGMA (auto-read from run_stats if empty)

# === Genotype input format ===
# Set to 'plink' to convert PLINK binary files to VCF automatically.
GENO_FORMAT="vcf"                 # vcf | plink
PLINK_PREFIX=""                   # e.g. /path/to/genotypes (without .bed/.bim/.fam)

# === Inputs ===
# VCF path (used directly when GENO_FORMAT=vcf; auto-set when GENO_FORMAT=plink)
# If the .tbi index is absent it is created automatically with tabix.
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
SLURM_EXTRA_ARGS=""           # extra "#SBATCH ..." line(s) (e.g. #SBATCH --account=...)

# === SLURM — merge job (submitted with afterok dependency) ===
SLURM_MERGE_TIME="02:00:00"   # wall-time for the merge step
SLURM_MERGE_MEM="16G"         # memory for the merge step
SLURM_MERGE_CPUS=2            # CPUs for the merge step

# === SLURM — steps IV-VI job (after all merge jobs) ===
SLURM_POST_TIME="04:00:00"
SLURM_POST_MEM="16G"
SLURM_POST_CPUS=4

# === Commands run at the start of every submitted job (module loads, PATH) ===
SLURM_SETUP=""
```

The chromosome job array is built automatically from the `CHROMOSOMES` variable, so setting e.g. `CHROMOSOMES="1 2 3"` will submit `--array=1,2,3` rather than the full 1–22.

Generated job scripts are self-contained (paths and settings are written into them, and
`SLURM_SETUP` runs first), because on some clusters — e.g. Snellius — jobs submitted from inside
a job do not inherit its environment.

---

## Running on Snellius (SURF) / HPC without internet on compute nodes

`examples/snellius/` contains a launcher, an example config and a Nextflow override that were used
to reproduce the ALFA regional-WMH sex-stratified analysis (merged results byte-identical to the
original run) and to run its AMYPAD replication:

```bash
cp examples/snellius/snellius.conf.example my.conf      # edit the <...> paths
bash   examples/snellius/run_snellius.sh check my.conf  # login node: steps I-II
sbatch examples/snellius/run_snellius.sh run   my.conf --skip-test
# after step III, on the login node (internet access), for the rsIDs of step IV:
bash scripts/make_dbsnp_subset.sh /path/dbsnp_top1000.vcf.gz 1000 OUT/results{,_female,_male}/mvgwas_merged.tsv
```

Points that matter there:
- **Home quota is 1M files** — put `OUTPUT_DIR` on project space (Nextflow writes ~10 files per chunk).
- **No internet on compute nodes** — set `DBSNP_VCF` to a local file (see `scripts/make_dbsnp_subset.sh`;
  remote region queries from Snellius are ~1 KB/s, a full download of dbSNP ~8 min).
- **Modules** — put them in `SLURM_SETUP`; the launcher and every submitted job run it.
- **Billing** — genoa jobs are charged for at least 24 CPUs.

### Verifying a reproduction

Compare the merged results with checksums from an earlier run (one `<md5> <stratum>` line per stratum):

```bash
bash examples/snellius/run_snellius.sh compare my.conf expected_md5.txt
# IDENTICAL  combined: OUT/results/mvgwas_merged.tsv
# IDENTICAL  male:     OUT/results_male/mvgwas_merged.tsv
# IDENTICAL  female:   OUT/results_female/mvgwas_merged.tsv
```

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

To re-run only post-processing, enrichment and the report:

```bash
bash run.sh --config my_run.conf --steps 4,5,6 --resume
```

With `--resume`, the run metadata and `run_stats.tsv` from earlier steps are kept (appended to), so
later steps still know the strata and per-stratum inputs.

---

## Input File Formats

### Phenotype / Covariate TSV
- Tab-separated, **header row required**
- First column must be **sample ID** (matching VCF sample names)
- All phenotype and covariate columns are passed directly to `mvgwas-nf`

### VCF
- Bgzipped + tabix-indexed (`.vcf.gz` + `.vcf.gz.tbi`)
- If the `.tbi` index is absent, it is created automatically using `tabix` — no manual step needed
- Sample IDs must overlap with phenotype / covariate files

### PLINK binary fileset (optional)
- Set `GENO_FORMAT=plink` and `PLINK_PREFIX=/path/to/prefix` in your config
- The pipeline expects `.bed`, `.bim`, and `.fam` files all sharing the same prefix
- Conversion to bgzipped VCF is performed automatically using `plink2` (preferred) or `plink 1.9` as a fallback, restricted to autosomes (chromosomes 1–22)
- The converted file is written to `OUTPUT_DIR/inputs/<prefix>.vcf.gz` and indexed; `VCF_FILE` is updated automatically for all downstream steps
- `plink2` or `plink` must be available in `PATH`

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Nextflow not found` | Install: `curl -s https://get.nextflow.io | bash` then move to PATH |
| `bcftools: command not found` | `conda install -c bioconda bcftools` or `brew install bcftools` |
| `plink2: command not found` | `conda install -c bioconda plink2` or `brew install plink2` (only needed for `GENO_FORMAT=plink`) |
| VCF index missing | Automatically recreated by the pipeline; or run `tabix -p vcf your.vcf.gz` manually |
| Java 8 detected | Auto-handled (NXF_VER=22.04.0); or upgrade Java |
| SLURM jobs queued forever | Check `SLURM_PARTITION` in config |
| rsID annotation slow, hangs, or all rsIDs empty | Set `DBSNP_VCF` to a local file (`scripts/make_dbsnp_subset.sh`); compute nodes often have no internet |
| `Disk quota exceeded` / `Unable to create file .command.run` | File-count quota reached: move `OUTPUT_DIR` to project/scratch space |
| `OUTPUT_DIR: unbound variable` in a SLURM job | Jobs do not inherit the environment on your cluster: put module loads in `SLURM_SETUP` (job scripts already contain all paths) |
| Steps IV–VI never ran in SLURM mode | A chromosome or merge job failed (`afterok` dependency): check `logs/chr_*`/`merge_*`, fix, then `bash run.sh --config ... --resume` |
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
