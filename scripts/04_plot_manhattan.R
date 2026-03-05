#!/usr/bin/env Rscript
# =============================================================================
# scripts/04_plot_manhattan.R  —  Manhattan, QQ, and Regional Manhattan plots
# =============================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(cowplot)
})

# =============================================================================
# CLI arguments
# =============================================================================
option_list <- list(
  make_option("--input",          type="character", help="Full merged GWAS TSV"),
  make_option("--top-file",       type="character", help="Top SNPs TSV (with optional rsID column)"),
  make_option("--output-dir",     type="character", help="Directory for output plots"),
  make_option("--run-name",       type="character", default="gwas_run", help="Run label for file names"),
  make_option("--pval-col",       type="character", default="P_manta", help="P-value column name"),
  make_option("--genome-build",   type="character", default="GRCh38",  help="GRCh37 or GRCh38"),
  make_option("--gw-thresh",      type="double",    default=5e-8,  help="Genome-wide significance threshold"),
  make_option("--sug-thresh",     type="double",    default=1e-5,  help="Suggestive threshold"),
  make_option("--top-label-n",    type="integer",   default=20,    help="Number of top SNPs to label"),
  make_option("--regional-window",type="integer",   default=2000000, help="Regional window (bp) around top locus"),
  make_option("--chr-col",        type="character", default=NULL,  help="Chromosome column name (auto-detected)"),
  make_option("--pos-col",        type="character", default=NULL,  help="Position column name (auto-detected)")
)
opts <- parse_args(OptionParser(option_list=option_list))

if (is.null(opts$input))      stop("--input is required")
if (is.null(opts[["top-file"]])) stop("--top-file is required")
if (is.null(opts[["output-dir"]])) stop("--output-dir is required")

dir.create(opts[["output-dir"]], recursive=TRUE, showWarnings=FALSE)
RUN   <- opts[["run-name"]]
PVCOL <- opts[["pval-col"]]
GW    <- opts[["gw-thresh"]]
SUG   <- opts[["sug-thresh"]]
NLBL  <- opts[["top-label-n"]]
WIN   <- opts[["regional-window"]]
BUILD <- opts[["genome-build"]]

cat(sprintf("[INFO] Loading merged results: %s\n", opts$input))
dt <- fread(opts$input, showProgress=FALSE)
cat(sprintf("[INFO] Loaded %d rows × %d cols\n", nrow(dt), ncol(dt)))

# ── Auto-detect column names ─────────────────────────────────────────────────
auto_col <- function(dt, candidates) {
  cols_lower <- tolower(names(dt))
  for (c in candidates) {
    idx <- which(cols_lower == tolower(c))
    if (length(idx)) return(names(dt)[idx[1]])
  }
  return(NULL)
}
CHR_COL <- opts[["chr-col"]] %||% auto_col(dt, c("CHR","#CHROM","CHROM","chromosome"))
POS_COL <- opts[["pos-col"]] %||% auto_col(dt, c("POS","BP","position"))

if (is.null(CHR_COL)) stop("Cannot find CHR column. Use --chr-col to specify.")
if (is.null(POS_COL)) stop("Cannot find POS column. Use --pos-col to specify.")
cat(sprintf("[INFO] CHR col='%s'  POS col='%s'  P col='%s'\n", CHR_COL, POS_COL, PVCOL))

if (!PVCOL %in% names(dt)) {
  cat(sprintf("[WARN] P-value column '%s' not found. Available: %s\n",
              PVCOL, paste(names(dt)[1:min(20,ncol(dt))], collapse=", ")))
  stop("P-value column not found.")
}

# ── Clean and prepare data ───────────────────────────────────────────────────
setnames(dt, c(CHR_COL, POS_COL, PVCOL), c("CHR","POS","P"))
dt[, CHR := as.integer(sub("^chr","",as.character(CHR)))]
dt[, POS := as.integer(POS)]
dt[, P   := as.numeric(P)]
dt <- dt[!is.na(CHR) & !is.na(POS) & !is.na(P) & P > 0 & P <= 1]
dt[, LOGP := -log10(P)]
cat(sprintf("[INFO] %d valid associations after cleaning\n", nrow(dt)))

# Assign cumulative genomic position for Manhattan
CHR_ORDER <- 1:22
dt <- dt[CHR %in% CHR_ORDER]
setorder(dt, CHR, POS)
chr_offsets <- dt[, .(maxpos = max(POS)), by=CHR]
setorder(chr_offsets, CHR)
chr_offsets[, offset := cumsum(as.numeric(c(0, maxpos[-.N]))) + seq_len(.N)*2e7]
dt <- merge(dt, chr_offsets[, .(CHR, offset)], by="CHR")
dt[, gpos := POS + offset]

chr_mids <- dt[, .(mid = mean(range(gpos))), by=CHR]

# ── Load top SNPs for labels ─────────────────────────────────────────────────
top_dt <- fread(opts[["top-file"]], showProgress=FALSE)
if (!PVCOL %in% names(top_dt)) {
  pcol2 <- auto_col(top_dt, c("P_manta","P","pvalue","p_value","Pvalue"))
  if (!is.null(pcol2)) setnames(top_dt, pcol2, PVCOL)
}
top_chrcol <- auto_col(top_dt, c("CHR","#CHROM","CHROM"))
top_poscol <- auto_col(top_dt, c("POS","BP","position"))
if (!is.null(top_chrcol)) setnames(top_dt, top_chrcol, "CHR")
if (!is.null(top_poscol)) setnames(top_dt, top_poscol, "POS")
top_dt[, CHR := as.integer(sub("^chr","",as.character(CHR)))]
top_dt[, POS := as.integer(POS)]
top_dt[, P   := as.numeric(get(PVCOL))]
top_dt[, LOGP:= -log10(P)]
top_dt <- merge(top_dt, chr_offsets[, .(CHR, offset)], by="CHR")
top_dt[, gpos := POS + offset]
top_dt <- head(top_dt[order(P)], NLBL)

# ── Build label column ───────────────────────────────────────────────────────
label_col <- NULL
for (lc in c("rsID","ID","SNP","rsid","snp","name")) {
  if (lc %in% names(top_dt)) { label_col <- lc; break }
}
if (!is.null(label_col)) {
  top_dt[, label := fifelse(
    is.na(get(label_col)) | get(label_col) %in% c("",".",NA_character_),
    paste0("chr", CHR, ":", POS),
    get(label_col)
  )]
} else {
  top_dt[, label := paste0("chr", CHR, ":", POS)]
}

# =============================================================================
# Plot 1: Manhattan
# =============================================================================
cat("[INFO] Generating Manhattan plot...\n")
dt[, chr_col := factor(CHR %% 2, labels=c("even","odd"))]

p_mht <- ggplot(dt, aes(x=gpos, y=LOGP, color=chr_col)) +
  geom_point(size=0.4, alpha=0.7, show.legend=FALSE) +
  scale_color_manual(values=c("even"="#4E79A7","odd"="#A0CBE8")) +
  geom_hline(yintercept=-log10(GW),  color="red",    linetype="dashed", linewidth=0.6) +
  geom_hline(yintercept=-log10(SUG), color="orange", linetype="dotted", linewidth=0.5) +
  scale_x_continuous(
    breaks = chr_mids$mid,
    labels = chr_mids$CHR,
    expand = expansion(mult=c(0.01,0.01))
  ) +
  scale_y_continuous(expand=expansion(mult=c(0.01,0.08))) +
  geom_point(data=top_dt, aes(x=gpos, y=LOGP), color="#E15759", size=1.5, inherit.aes=FALSE) +
  geom_label_repel(
    data    = top_dt,
    aes(x=gpos, y=LOGP, label=label),
    size            = 2.8,
    color           = "#222222",
    box.padding     = 0.4,
    point.padding   = 0.3,
    max.overlaps    = 25,
    segment.size    = 0.3,
    segment.alpha   = 0.6,
    min.segment.len = 0.2,
    inherit.aes     = FALSE
  ) +
  annotate("text", x=max(dt$gpos)*0.97, y=-log10(GW)+0.2,
           label=sprintf("GW: p=%.0e", GW), size=2.5, hjust=1, color="red") +
  labs(
    title   = sprintf("Manhattan Plot — %s", RUN),
    x       = "Chromosome",
    y       = expression(-log[10](italic(P)))
  ) +
  theme_classic(base_size=11) +
  theme(
    plot.title      = element_text(hjust=0.5, size=13, face="bold"),
    axis.text.x     = element_text(size=7),
    panel.grid.major.y = element_line(color="grey90", linewidth=0.3)
  )

ggsave(
  file.path(opts[["output-dir"]], sprintf("manhattan_%s.png", RUN)),
  p_mht, width=14, height=5, dpi=200
)
cat(sprintf("[OK] Manhattan saved\n"))

# =============================================================================
# Plot 2: QQ plot
# =============================================================================
cat("[INFO] Generating QQ plot...\n")

# Lambda GC
chisq <- qchisq(dt$P, df=1, lower.tail=FALSE)
lambda <- round(median(chisq, na.rm=TRUE) / qchisq(0.5, df=1), 3)

n  <- nrow(dt)
observed  <- sort(dt$LOGP, decreasing=TRUE)
expected  <- -log10(ppoints(n))
qq_df <- data.frame(expected=expected, observed=observed)

p_qq <- ggplot(qq_df, aes(x=expected, y=observed)) +
  geom_point(size=0.5, alpha=0.6, color="#4E79A7") +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  annotate("text", x=0.1, y=max(observed)*0.95,
           label=sprintf("\u03bb = %.3f", lambda), hjust=0, size=4) +
  labs(
    title = sprintf("QQ Plot — %s", RUN),
    x     = expression(Expected~-log[10](italic(P))),
    y     = expression(Observed~-log[10](italic(P)))
  ) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(hjust=0.5, face="bold"))

ggsave(
  file.path(opts[["output-dir"]], sprintf("qq_%s.png", RUN)),
  p_qq, width=6, height=6, dpi=200
)
cat(sprintf("[OK] QQ plot saved  (lambda=%.3f)\n", lambda))

# =============================================================================
# Plot 3: Regional Manhattan around top locus
# =============================================================================
cat("[INFO] Generating regional Manhattan plot...\n")

top1  <- top_dt[which.min(P)]
t_chr <- top1$CHR
t_pos <- top1$POS
t_lbl <- top1$label
HALF  <- WIN / 2

cat(sprintf("[INFO] Top locus: %s  chr%d:%d\n", t_lbl, t_chr, t_pos))

# Subset data to window
reg <- dt[CHR == t_chr & POS >= (t_pos - HALF) & POS <= (t_pos + HALF)]
cat(sprintf("[INFO] Regional: %d variants in window\n", nrow(reg)))

if (nrow(reg) < 5) {
  cat("[WARN] Too few variants in region; skipping regional plot\n")
} else {
  # ── Fetch gene annotations from UCSC REST API ────────────────────────────
  genes_df <- NULL
  tryCatch({
    ucsc_build <- if (BUILD == "GRCh38") "hg38" else "hg19"
    url <- sprintf(
      "https://api.genome.ucsc.edu/getData/track?genome=%s&track=ncbiRefSeq&chrom=chr%d&start=%d&end=%d",
      ucsc_build, t_chr, t_pos - HALF, t_pos + HALF
    )
    cat(sprintf("[INFO] Fetching gene track: %s\n", url))
    resp <- tryCatch(
      jsonlite::fromJSON(url, simplifyDataFrame=TRUE),
      error = function(e) NULL
    )
    if (!is.null(resp) && "ncbiRefSeq" %in% names(resp)) {
      g <- as.data.frame(resp$ncbiRefSeq)
      if (nrow(g) > 0) {
        g <- g[, intersect(c("name2","chromStart","chromEnd","strand"), names(g))]
        g <- g[!duplicated(g$name2), ]
        genes_df <- g
        cat(sprintf("[INFO] Got %d gene annotations\n", nrow(genes_df)))
      }
    }
  }, error = function(e) {
    cat(sprintf("[WARN] Could not fetch gene track: %s\n", conditionMessage(e)))
  })

  # ── Build regional plot (scatter panel) ──────────────────────────────────
  reg_label <- reg[which.min(P)]
  reg_label[, label := t_lbl]

  p_reg_scatter <- ggplot(reg, aes(x=POS/1e6, y=LOGP)) +
    geom_point(size=1.2, alpha=0.8, color="#4E79A7") +
    geom_point(data=reg_label, color="red", size=3, inherit.aes=FALSE,
               aes(x=POS/1e6, y=LOGP)) +
    geom_label_repel(data=reg_label, aes(x=POS/1e6, y=LOGP, label=label),
                     color="red", size=3, box.padding=0.5, inherit.aes=FALSE) +
    geom_hline(yintercept=-log10(GW),  color="red",    linetype="dashed", linewidth=0.5) +
    geom_hline(yintercept=-log10(SUG), color="orange", linetype="dotted", linewidth=0.4) +
    scale_x_continuous(labels=comma_format()) +
    labs(
      title = sprintf("Regional Manhattan — %s  (chr%d: %.2f–%.2f Mb)",
                      t_lbl, t_chr, (t_pos-HALF)/1e6, (t_pos+HALF)/1e6),
      x     = sprintf("chr%d position (Mb)", t_chr),
      y     = expression(-log[10](italic(P)))
    ) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(hjust=0.5, face="bold"))

  # ── Gene track panel (if available) ──────────────────────────────────────
  if (!is.null(genes_df) && nrow(genes_df) > 0) {
    names(genes_df)[names(genes_df) == "chromStart"] <- "start"
    names(genes_df)[names(genes_df) == "chromEnd"]   <- "end"
    names(genes_df)[names(genes_df) == "name2"]      <- "gene"
    genes_df$y_level <- as.numeric(factor(genes_df$gene)) %% 5

    p_genes <- ggplot(genes_df) +
      geom_segment(aes(x=start/1e6, xend=end/1e6, y=y_level, yend=y_level),
                   linewidth=2.5, color="#555555") +
      geom_text(aes(x=(start+end)/2e6, y=y_level+0.3, label=gene),
                size=2.5, hjust=0.5) +
      scale_x_continuous(
        limits = c((t_pos-HALF)/1e6, (t_pos+HALF)/1e6),
        labels = comma_format()
      ) +
      labs(x=sprintf("chr%d (Mb)", t_chr), y="") +
      theme_classic(base_size=9) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            axis.line.y=element_blank())

    p_reg <- plot_grid(p_reg_scatter, p_genes, ncol=1,
                       rel_heights=c(3,1), align="v", axis="lr")
  } else {
    p_reg <- p_reg_scatter
  }

  ggsave(
    file.path(opts[["output-dir"]], sprintf("regional_%s.png", RUN)),
    p_reg, width=10, height=6, dpi=200
  )
  cat(sprintf("[OK] Regional Manhattan saved\n"))
}

cat(sprintf("[INFO] All plots written to: %s\n", opts[["output-dir"]]))
