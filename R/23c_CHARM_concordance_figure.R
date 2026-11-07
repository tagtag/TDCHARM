# ============================================================
# 23c_CHARM_concordance_figure.R
# GSE303006 / GENCODE M23
#
# Revision figure (Reviewer 1, Comment 5): concordance of core #12
# with the cell-type-resolved E-P links of the original CHARM study,
# drawn from the outputs of 23b (no new computation).
#
#   a  Core-#12 score (percentile among all 730,969 candidate pairs) of
#      candidate pairs that CHARM links only in the group-A cell type
#      (Ndnf/Lamp5), only in the group-B cell type (L3/4 IT), or only
#      in other cell types; Wilcoxon test A-only vs B-only (from 23b).
#   b  Number of group-A and group-B CHARM links among the 1,000
#      highest-scoring pairs on the positive and on the negative side
#      of the core-#12 ranking, with the number expected from their
#      overall frequency; common odds ratio and P of the CMH test
#      stratified by enhancer-TSS distance decile (from 23b).
#
# Inputs (working directory): outputs of 23b, metadata
# Outputs:
#   GSE303006_M23_rev23c_CHARM_concordance.pdf / .png
#   GSE303006_M23_rev23c_log.txt
# ============================================================

library(data.table)

local({
    lib <- file.path(Sys.getenv("TDCHARM_REPO", unset = "."), "R", "revision_lib.R")
    if (!exists("rev_find_repo")) {
        if (!file.exists(lib)) {
            fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
            if (length(fa)) lib <- file.path(dirname(normalizePath(sub("^--file=", "", fa[1]))), "revision_lib.R")
            if (!file.exists(lib) && exists("repo", envir = globalenv()))
                lib <- file.path(get("repo", envir = globalenv()), "R", "revision_lib.R")
        }
        sys.source(lib, envir = globalenv())
    }
})

TOP_N <- as.integer(strsplit(Sys.getenv("TDCHARM_TOP_N", unset = "1000,5000"), ",")[[1]])[1]
UNIVERSE <- "all_pairs"
SEED <- 20261105L

PAIRS_FILE <- "GSE303006_M23_rev23b_CHARM_linked_pairs.tsv.gz"
AB_FILE <- "GSE303006_M23_rev23b_groupA_vs_groupB_links.tsv"
CMH_FILE <- "GSE303006_M23_rev23b_CMH_topN.tsv"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_PDF <- "GSE303006_M23_rev23c_CHARM_concordance.pdf"
OUT_PNG <- "GSE303006_M23_rev23c_CHARM_concordance.png"
OUT_LOG <- "GSE303006_M23_rev23c_log.txt"

need <- c(PAIRS_FILE, AB_FILE, CMH_FILE, META_FILE)
if (any(!file.exists(need))) stop("Missing (run 23b first): ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "23c_CHARM_concordance_figure.R")

pairs <- fread(PAIRS_FILE)
ab <- fread(AB_FILE)
cmh <- fread(CMH_FILE)
groups <- rev_contrast_groups(fread(META_FILE)$cell_type)
short <- function(x) gsub("_", " ", sub("^(Inh|Ex)[_ ]+", "", x))
labA <- paste(short(groups$A), collapse = "/")
labB <- paste(short(groups$B), collapse = "/")
cat("Group A: ", labA, "; group B: ", labB, "\n", sep = "")

COL_A <- "#B2182B"   # positive side (red, as in Fig. 4)
COL_B <- "#2166AC"   # negative side (blue)
COL_O <- "grey55"

fmt_p <- function(p) {
    if (!is.finite(p)) return("NA")
    if (p >= 1e-3) return(sprintf("P = %.3f", p))
    e <- floor(log10(p)); m <- p / 10^e
    bquote(P == .(sprintf("%.1f", m)) %*% 10^.(e))
}

# ------------------------------------------------------------
# Panel data
# ------------------------------------------------------------

cls_order <- c("A_only", "B_only", "other")
cls_lab <- c(paste0(labA, "\nonly"), paste0(labB, "\nonly"), "other cell\ntypes only")
cls_col <- c(COL_A, COL_B, COL_O)
pa <- pairs[cls %in% cls_order]
pa[, cls := factor(cls, cls_order)]
cat("\nPanel a: CHARM-linked candidate pairs by class\n")
print(pa[, .(n = .N, median_percentile = median(s_percentile),
             frac_positive = mean(s > 0)), by = cls][order(cls)])
print(ab)

pb <- cmh[universe == UNIVERSE & top_n == TOP_N &
              reference %in% c("CHARM_groupA", "CHARM_groupB") &
              ranking %in% c("positive", "negative")]
pb[, expected := fraction_overall * top_n]
pb[, reference := factor(reference, c("CHARM_groupA", "CHARM_groupB"))]
pb[, ranking := factor(ranking, c("positive", "negative"))]
setorder(pb, reference, ranking)
cat("\nPanel b: top-", TOP_N, " pairs (universe ", UNIVERSE, ")\n", sep = "")
print(pb[, .(reference, ranking, supported_in_top, expected = round(expected, 2),
             OR = round(common_OR, 2), CMH_p = signif(CMH_p, 3))])
has_b <- nrow(pb) == 4L
if (!has_b) warning("Top-", TOP_N, " CMH results not available in 23b output; panel b left empty")

# ------------------------------------------------------------
# Drawing
# ------------------------------------------------------------

draw <- function() {
    set.seed(SEED)
    layout(matrix(1:2, 1), widths = c(1.1, 1))

    # a
    par(mar = c(5.2, 4.6, 3.2, 1), mgp = c(2.6, 0.7, 0), las = 1)
    plot(NA, xlim = c(0.5, 3.5), ylim = c(0, 108), xaxt = "n", yaxt = "n",
         xlab = "", ylab = "Core #12 score percentile\n(all candidate E-P pairs)", bty = "l")
    axis(2, at = seq(0, 100, 25))
    abline(h = 50, lty = 2, col = "grey70")
    for (i in seq_along(cls_order)) {
        y <- pa[cls == cls_order[i], s_percentile]
        if (!length(y)) next
        points(i + runif(length(y), -0.28, 0.28), y, pch = 16, cex = 0.45,
               col = adjustcolor(cls_col[i], 0.45))
        boxplot(y, at = i, add = TRUE, boxwex = 0.55, outline = FALSE, axes = FALSE,
                col = NA, border = "black", lwd = 1.1)
    }
    axis(1, at = 1:3, labels = FALSE)
    n_cls <- pa[, .N, keyby = cls][J(cls_order), N]
    n_cls[is.na(n_cls)] <- 0L
    text(1:3, par("usr")[3] - 4, paste0(cls_lab, "\n(n = ", n_cls, ")"),
         xpd = NA, adj = c(0.5, 1), cex = 0.85)
    if ("rank_biserial_A_minus_B" %in% names(ab)) {
        segments(1, 103, 2, 103); segments(c(1, 2), 103, c(1, 2), 100)
        text(1.5, 106.5, bquote(r[rb] == .(sprintf("%.2f", ab$rank_biserial_A_minus_B)) * ", " *
                                    .(fmt_p(ab$wilcox_p))), cex = 0.8, xpd = NA)
    }
    mtext("a", side = 3, line = 1.6, adj = -0.25, font = 2, cex = 1.2)

    # b
    par(mar = c(5.2, 4.6, 3.2, 1))
    if (!has_b) {
        plot.new(); text(0.5, 0.5, paste0("top-", TOP_N, " results\nnot available"), cex = 0.9)
        mtext("b", side = 3, line = 1.6, adj = -0.25, font = 2, cex = 1.2)
        return(invisible())
    }
    m <- matrix(pb$supported_in_top, 2, dimnames = list(c("positive", "negative"), c("A", "B")))
    ymax <- max(m, 1) * 1.75
    bp <- barplot(m, beside = TRUE, col = c(COL_A, COL_B), border = NA, ylim = c(0, ymax),
                  space = c(0.15, 1.3), names.arg = rep("", 2), xlim = c(0.3, 7.4),
                  ylab = paste0("CHARM links among top ", format(TOP_N, big.mark = ","), " pairs"))
    text(colMeans(bp), par("usr")[3] - ymax * 0.04,
         c(paste0(labA, "\nlinks"), paste0(labB, "\nlinks")), xpd = NA, adj = c(0.5, 1), cex = 0.85)
    e <- matrix(pb$expected, 2)
    segments(bp - 0.45, e, bp + 0.45, e, lwd = 2)
    for (j in 1:4) {
        lab <- if (pb$supported_in_top[j] == 0) "none" else
            bquote(atop("OR" == .(sprintf("%.0f", pb$common_OR[j])), .(fmt_p(pb$CMH_p[j]))))
        text(bp[j], max(m[j], e[j]) + ymax * 0.03, lab, adj = c(0.5, 0), cex = 0.62, xpd = NA)
    }
    legend("topleft", inset = c(0.02, 0), bty = "n", cex = 0.72,
           fill = c(COL_A, COL_B, NA), border = c(NA, NA, NA),
           legend = c("positive side", "negative side", "expected"),
           lty = c(NA, NA, 1), lwd = c(NA, NA, 2), merge = TRUE)
    mtext("b", side = 3, line = 1.6, adj = -0.25, font = 2, cex = 1.2)
}

W <- 7.4; H <- 4.2
pdf(OUT_PDF, width = W, height = H, useDingbats = FALSE); draw(); invisible(dev.off())
png(OUT_PNG, width = W * 300, height = H * 300, res = 300); draw(); invisible(dev.off())
cat("\nSaved ", OUT_PDF, " and ", OUT_PNG, "\n", sep = "")
rev_end_log(log_con)
