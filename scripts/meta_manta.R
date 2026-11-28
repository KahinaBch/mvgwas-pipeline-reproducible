#!/usr/bin/env Rscript
## =============================================================================
## meta_manta.R — fixed-effect multivariate (MANTA) meta-analysis across cohorts
## -----------------------------------------------------------------------------
## Input per cohort (from 04b_meta_export.R): <cohort>.mvsumstats.tsv.gz and
## <cohort>.rescor.tsv. For each variant:
##   1. align BETA to a common effect allele (sign flip if REF/ALT swapped);
##   2. inverse-variance meta-analysis of each trait's BETA;
##   3. MANTA statistic on the pooled betas, T = sum_j beta_j^2, whose null
##      distribution is sum_i lambda_i chi2_1 with lambda = eig(D R D),
##      D = diag(pooled SE), R = pooled residual trait correlation.
## Also reports the Wald / MANOVA-type statistic z' R^-1 z ~ chi2_q (P_wald_meta).
## For a single cohort this reproduces that cohort's P_manta (up to R being
## estimated without the SNP in the model).
##
## Usage:
##   Rscript meta_manta.R --sumstats A.mvsumstats.tsv.gz,B.mvsumstats.tsv.gz \
##                        --rescor   A.rescor.tsv,B.rescor.tsv \
##                        --out meta_manta.tsv.gz [--min-cohorts 2 (default)]
## =============================================================================
suppressPackageStartupMessages({ library(optparse); library(data.table) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sumstats",    type = "character", help = "comma-separated *.mvsumstats.tsv.gz"),
  make_option("--rescor",      type = "character", help = "comma-separated *.rescor.tsv (same order)"),
  make_option("--out",         type = "character", help = "output file"),
  make_option("--min-cohorts", type = "integer", default = 2L, help = "min cohorts per variant [%default]")
)))
ss_f <- strsplit(opt$sumstats, ",")[[1]]; rc_f <- strsplit(opt$rescor, ",")[[1]]
stopifnot(length(ss_f) == length(rc_f))
K <- length(ss_f)

## p-value of sum_i lambda_i chi2_1 > T, by Farebrother's AS 204 algorithm as in
## manta::p.asympt (eps starts at 1e-14, relaxed x2 until it converges, p floored
## at eps). Uses manta when installed, else CompQuadForm::farebrother (same algorithm).
pfun <- if (requireNamespace("manta", quietly = TRUE)) {
  function(T, lam) manta:::p.asympt(ss = T, df = 1, lambda = lam)[1]
} else if (requireNamespace("CompQuadForm", quietly = TRUE)) {
  function(T, lam, eps = 1e-14) {
    repeat {
      r <- CompQuadForm::farebrother(T, lam, eps = eps)
      if (r$ifault %in% c(0, 4) || eps > 1e-10) break
      eps <- eps * 2
    }
    max(r$Qq, eps)
  }
} else stop("Install 'CompQuadForm' (install.packages('CompQuadForm')) or 'manta'", call. = FALSE)

ss <- lapply(ss_f, fread)
tr_k   <- lapply(ss, function(d) sub("^BETA_", "", grep("^BETA_", names(d), value = TRUE)))
traits <- Reduce(intersect, tr_k)
dropped <- setdiff(Reduce(union, tr_k), traits)
if (length(dropped)) warning("Traits not present in every cohort, left out: ", paste(dropped, collapse = ", "), call. = FALSE)
if (!length(traits)) stop("No trait is shared by all cohorts", call. = FALSE)
q <- length(traits)
cat(sprintf("[INFO] %d cohorts, %d shared traits: %s\n", K, q, paste(traits, collapse = ", ")))

## Pooled residual correlation, weighted by cohort N (median N over variants)
Rs <- lapply(rc_f, function(f) { r <- data.frame(fread(f), row.names = 1, check.names = FALSE)
                                  as.matrix(r[traits, traits]) })
Ns <- vapply(ss, function(d) median(d$N, na.rm = TRUE), 1)
R  <- Reduce(`+`, Map(`*`, Rs, Ns)) / sum(Ns)

## Align: effect allele = MARKER's 4th field (alphabetically larger allele)
markers <- unique(unlist(lapply(ss, `[[`, "MARKER")))
idx <- lapply(ss, function(d) match(markers, d$MARKER))
SW <- SWB <- matrix(0, length(markers), q); ncoh <- integer(length(markers)); Ntot <- numeric(length(markers))
for (k in seq_len(K)) {
  d <- ss[[k]][idx[[k]]]
  ea <- sub(".*:", "", markers)
  sgn <- ifelse(d$ALT == ea, 1, -1)
  B <- as.matrix(d[, paste0("BETA_", traits), with = FALSE]) * sgn
  S <- as.matrix(d[, paste0("SE_",   traits), with = FALSE])
  W <- 1 / S^2; ok <- is.finite(B) & is.finite(W) & S > 0
  W[!ok] <- 0; B[!ok] <- 0
  SW <- SW + W; SWB <- SWB + W * B
  has <- rowSums(ok) == q; ncoh <- ncoh + has; Ntot <- Ntot + ifelse(has, d$N, 0)
}
keep <- ncoh >= opt[["min-cohorts"]] & rowSums(SW > 0) == q
bm <- (SWB / SW)[keep, , drop = FALSE]; sm <- sqrt(1 / SW)[keep, , drop = FALSE]

Rinv <- solve(R)
PP <- vapply(seq_len(nrow(bm)), function(i) {
  C   <- R * tcrossprod(sm[i, ])                          # Cov(pooled betas) = D R D
  lam <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  lam <- lam[lam / sum(lam) > 1e-3]                       # same tolerance as manta
  z   <- bm[i, ] / sm[i, ]
  c(pfun(sum(bm[i, ]^2), lam),                            # MANTA: beta' beta
    pchisq(drop(crossprod(z, Rinv %*% z)), q, lower.tail = FALSE))  # Wald: beta' C^-1 beta
}, c(1, 1))
P <- PP[1, ]; Pw <- PP[2, ]

f <- tstrsplit(markers[keep], ":", fixed = TRUE)
res <- data.table(MARKER = markers[keep], CHR = f[[1]], POS = as.integer(f[[2]]),
                  EA = f[[4]], NEA = f[[3]], N_cohorts = ncoh[keep], N = Ntot[keep],
                  P_manta_meta = P, P_wald_meta = Pw)
res[, paste0("BETA_", traits) := as.data.table(bm)]
res[, paste0("SE_",   traits) := as.data.table(sm)]
res[, paste0("P_",    traits) := as.data.table(2 * pnorm(-abs(bm / sm)))]
res <- res[order(suppressWarnings(as.integer(CHR)), CHR, POS)]
fwrite(res, opt$out, sep = "\t")
cat(sprintf("[INFO] %d variants written to %s\n", nrow(res), opt$out))
