#!/usr/bin/env Rscript
## =============================================================================
## 04b_meta_export.R — export one cohort's results for cross-cohort meta-analysis
## -----------------------------------------------------------------------------
## Needs a run with PER_TRAIT_STATS=true (mvgwas-nf --pt), so that
## mvgwas_merged.tsv carries N, AF_ALT, BETA_<trait>, SE_<trait>.
##
## Writes to --out-dir:
##   <cohort>.rescor.tsv          trait x trait residual correlation (Y ~ covariates)
##   <cohort>.mvsumstats.tsv.gz   compact wide file (MARKER, alleles, N, BETA_*, SE_*)
##                                used by meta_manta.R for the multivariate meta-analysis
## and, only with --metal true (META_WRITE_METAL=true; ~120 bytes x variants x traits):
##   <cohort>.<trait>.metal.tsv   one METAL-ready file per trait
##   metal_<trait>.txt            METAL script stub for that trait (add other cohorts)
## =============================================================================
suppressPackageStartupMessages({ library(optparse); library(data.table) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--merged",    type = "character", help = "mvgwas_merged.tsv"),
  make_option("--pheno",     type = "character", help = "phenotype file used for the run"),
  make_option("--cov",       type = "character", help = "covariate file used for the run"),
  make_option("--transform", type = "character", default = "none", help = "none | sqrt | log [%default]"),
  make_option("--cohort",    type = "character", help = "cohort / stratum label used in file names"),
  make_option("--out-dir",   type = "character", help = "output directory"),
  make_option("--metal",     type = "logical", default = FALSE,
              help = "also write one METAL input file + script stub per trait [%default]")
)))
for (k in c("merged", "pheno", "cov", "cohort", "out-dir"))
  if (is.null(opt[[k]])) stop("missing --", k, call. = FALSE)
dir.create(opt[["out-dir"]], recursive = TRUE, showWarnings = FALSE)
od <- function(f) file.path(opt[["out-dir"]], f)

## ---- 1. per-trait files ------------------------------------------------------
dt <- fread(opt$merged)
bcols <- grep("^BETA_", names(dt), value = TRUE)
if (!length(bcols) || !all(c("N", "AF_ALT") %in% names(dt)))
  stop("No per-trait columns in ", opt$merged, " - rerun Step III with PER_TRAIT_STATS=true", call. = FALSE)
traits <- sub("^BETA_", "", bcols)

## Marker name independent of REF/ALT order and chr prefix, so the same variant
## gets the same name in every cohort; METAL then aligns alleles itself.
chr <- sub("^chr", "", as.character(dt$CHR))
a1 <- pmin(dt$REF, dt$ALT); a2 <- pmax(dt$REF, dt$ALT)
dt[, MARKER := paste(chr, POS, a1, a2, sep = ":")]

## METAL files are optional (by-hand heterogeneity checks): they can be large,
## and everything in them is also in the mvsumstats file below.
if (isTRUE(opt$metal)) for (tr in traits) {
  b <- dt[[paste0("BETA_", tr)]]; s <- dt[[paste0("SE_", tr)]]
  out <- data.table(MARKER = dt$MARKER, CHR = chr, POS = dt$POS,
                    EA = dt$ALT, NEA = dt$REF,          # BETA is per ALT allele
                    BETA = b, SE = s, P = 2 * pnorm(-abs(b / s)),
                    N = dt$N, EAF = dt$AF_ALT)
  fwrite(out[is.finite(BETA) & is.finite(SE) & SE > 0], od(sprintf("%s.%s.metal.tsv", opt$cohort, tr)), sep = "\t")
  writeLines(c(
    "SCHEME STDERR", "MARKER MARKER", "ALLELE EA NEA", "EFFECT BETA", "STDERR SE",
    "PVALUE P", "WEIGHT N", "FREQ EAF", "AVERAGEFREQ ON", "MINMAXFREQ ON",
    "SEPARATOR TAB", "",
    sprintf("PROCESS %s.%s.metal.tsv", opt$cohort, tr),
    "# PROCESS <other_cohort>.<trait>.metal.tsv", "",
    sprintf("OUTFILE META_%s_ .tbl", tr),
    "ANALYZE HETEROGENEITY", "QUIT"), od(sprintf("metal_%s.txt", tr)))
}

fwrite(dt[, c("MARKER", "CHR", "POS", "REF", "ALT", "N", "AF_ALT",
              paste0("BETA_", traits), paste0("SE_", traits)), with = FALSE],
       od(sprintf("%s.mvsumstats.tsv.gz", opt$cohort)), sep = "\t")

## ---- 2. residual trait correlation ------------------------------------------
## Same model as the GWAS without the SNP: transformed Y ~ covariates.
ph <- data.frame(fread(opt$pheno), row.names = 1, check.names = FALSE)
cv <- data.frame(fread(opt$cov),   row.names = 1, check.names = FALSE, stringsAsFactors = TRUE)
ids <- intersect(rownames(ph), rownames(cv))
Y <- as.matrix(ph[ids, traits, drop = FALSE])
Y <- switch(opt$transform, none = Y, sqrt = sqrt(Y), log = log(Y))
cv <- cv[ids, vapply(cv, function(x) length(unique(na.omit(x))) > 1, TRUE), drop = FALSE]  # as preprocess.R
cc <- complete.cases(Y, cv)
R <- cor(qr.resid(qr(model.matrix(~ ., cv[cc, , drop = FALSE])), Y[cc, , drop = FALSE]))
fwrite(data.table(trait = traits, R), od(sprintf("%s.rescor.tsv", opt$cohort)), sep = "\t")

cat(sprintf("[INFO] %s: %d variants, %d traits, N(rescor) = %d, METAL files: %s -> %s\n",
            opt$cohort, nrow(dt), length(traits), sum(cc),
            if (isTRUE(opt$metal)) "yes" else "no", opt[["out-dir"]]))
