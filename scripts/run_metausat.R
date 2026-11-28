#!/usr/bin/env Rscript
## =============================================================================
## run_metausat.R — metaUSAT / metaMANOVA (Ray & Boehnke 2018) on the per-cohort
## exports of 04b_meta_export.R (<cohort>.mvsumstats.tsv.gz, <cohort>.rescor.tsv)
## -----------------------------------------------------------------------------
## --mode pooled  (default): IVW per trait across cohorts -> K z-scores (assumes
##                homogeneous effects, independent cohorts; = README note 5)
## --mode stacked : all J*K cohort-trait z-scores jointly (allows heterogeneity)
## --rmode resid  (default): R from covariate residuals (block-diagonal if stacked)
## --rmode null   : R from z-scores of SNPs with all p > 1e-5 (handles overlap;
##                  needs many independent SNPs, i.e. genome-wide input)
## Fixes two issues in metaUSAT_v1.17.R: vector `weights` errors in R >= 4.2,
## and cor.pearson() returns NA when no SNP passes the p threshold.
## Needs: metaUSAT_v1.17.R (github.com/RayDebashree/metaUSAT), CompQuadForm,
##        minqa, psych, survey, data.table, optparse
## =============================================================================
suppressPackageStartupMessages({ library(optparse); library(data.table) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstats", type = "character", help = "comma-separated *.mvsumstats.tsv.gz"),
  make_option("--rescor",   type = "character", help = "comma-separated *.rescor.tsv (same order)"),
  make_option("--metausat", type = "character", help = "path to metaUSAT_v1.17.R"),
  make_option("--mode",     type = "character", default = "pooled"),
  make_option("--rmode",    type = "character", default = "resid"),
  make_option("--out",      type = "character")
)))
source(opt$metausat)
if (!exists("pchisqsum")) suppressPackageStartupMessages(library(survey))
i_w <- which(vapply(as.list(body(metausat)), function(e) any(grepl("weights != 1", deparse(e), fixed = TRUE)), TRUE))
stopifnot(length(i_w) == 1)
body(metausat)[[i_w]] <- quote(if (any(weights != 1)) {                   # fix 1
  if (length(weights) != k) stop("Length of vector 'weights' must be same as length of Z-vector.")
  Z <- Z * weights; W <- diag(weights); R <- W %*% R %*% W })

ss_f <- strsplit(opt$sumstats, ",")[[1]]; rc_f <- strsplit(opt$rescor, ",")[[1]]; J <- length(ss_f)
ss <- lapply(ss_f, fread)
traits <- Reduce(intersect, lapply(ss, function(d) sub("^BETA_", "", grep("^BETA_", names(d), value = TRUE))))
K <- length(traits)
Rs <- lapply(rc_f, function(f) as.matrix(data.frame(fread(f), row.names = 1, check.names = FALSE)[traits, traits]))
mk <- Reduce(intersect, lapply(ss, `[[`, "MARKER"))                  # variants in all cohorts
ea <- sub(".*:", "", mk)
Z <- SE <- B <- vector("list", J)
for (j in 1:J) {
  d <- ss[[j]][match(mk, MARKER)]; s <- ifelse(d$ALT == ea, 1, -1)   # align to common EA
  B[[j]]  <- s * as.matrix(d[, paste0("BETA_", traits), with = FALSE])
  SE[[j]] <- as.matrix(d[, paste0("SE_", traits), with = FALSE])
  Z[[j]]  <- B[[j]] / SE[[j]]
}
if (opt$mode == "pooled") {
  W  <- Reduce(`+`, lapply(SE, function(s) 1 / s^2))
  Zm <- Reduce(`+`, Map(function(b, s) b / s^2, B, SE)) / sqrt(W)
  N  <- vapply(ss, function(d) median(d$N, na.rm = TRUE), 1)
  R  <- Reduce(`+`, Map(`*`, Rs, N)) / sum(N)
} else {
  Zm <- do.call(cbind, Z)
  R  <- as.matrix(Matrix::bdiag(Rs))
}
if (opt$rmode == "null") {                                            # fix 2
  keep <- apply(2 * pnorm(-abs(Zm)), 1, function(p) all(p >= 1e-5))
  R <- cor(Zm[keep, , drop = FALSE])
}
res <- t(vapply(seq_len(nrow(Zm)), function(i) {
  o <- metausat(Zm[i, ], R, metamanova = TRUE)
  c(o$p.metamanova, o$p.metausat, o$omg.opt)
}, numeric(3)))
f <- tstrsplit(mk, ":", fixed = TRUE)
out <- data.table(MARKER = mk, CHR = f[[1]], POS = as.integer(f[[2]]), EA = f[[4]], NEA = f[[3]],
                  P_metaMANOVA = res[, 1], P_metaUSAT = res[, 2], omega_opt = res[, 3])
fwrite(out[order(suppressWarnings(as.integer(CHR)), CHR, POS)], opt$out, sep = "\t")
cat(sprintf("[INFO] %s mode, R from %s: %d variants, %d cohorts x %d traits -> %s\n",
            opt$mode, opt$rmode, nrow(out), J, K, opt$out))
