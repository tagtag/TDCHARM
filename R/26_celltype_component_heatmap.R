# ============================================================
# 26_celltype_component_heatmap.R
# GSE303006 / GENCODE M23
#
# Revision figure (Reviewer 2, Comment 2).
#
# Characterizes ALL cell components that occur among the 100
# largest core elements (not only CELL03) in one figure:
#
#   rows    : cell types (ordered by hierarchical clustering of the
#             profile; cell count shown at the right)
#   columns : cell components in the top-100 core elements
#             (number of top-100 core elements and cell-type
#             epsilon^2 shown above each column)
#   colour  : mean standardized loading of the cells of that type
#             (loading z-scored over all 4,258 cells per component),
#             diverging palette, clipped at +/- CLIP for display only
#
# The sign of each component is arbitrary (as for any singular
# vector); the figure shows which cell types are contrasted.
#
# Inputs (working directory): outputs of 08, 09, 10
# Outputs:
#   GSE303006_M23_rev26_celltype_by_component_mean_z.tsv
#   GSE303006_M23_rev26_celltype_by_component_heatmap.pdf
#   GSE303006_M23_rev26_celltype_by_component_heatmap.png
#   GSE303006_M23_rev26_log.txt
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

TOP_CORE <- 100L
CLIP <- 2
MIN_CELLS <- 10L
HIGHLIGHT <- "CELL03"

SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_TSV <- "GSE303006_M23_rev26_celltype_by_component_mean_z.tsv"
OUT_PDF <- "GSE303006_M23_rev26_celltype_by_component_heatmap.pdf"
OUT_PNG <- "GSE303006_M23_rev26_celltype_by_component_heatmap.png"
OUT_LOG <- "GSE303006_M23_rev26_log.txt"

need <- c(SMALL_FILE, TOP_CORE_FILE, META_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "26_celltype_component_heatmap.R")

SF <- readRDS(SMALL_FILE)
U <- SF$U_CELL
meta <- fread(META_FILE)[match(rownames(U), cell)]
top <- fread(TOP_CORE_FILE)[core_rank <= TOP_CORE]
nb <- top[, .N, by = b_CELL][order(b_CELL)]
comps <- sprintf("CELL%02d", nb$b_CELL)

Z <- scale(U[, nb$b_CELL, drop = FALSE])
colnames(Z) <- comps
ct <- meta$cell_type
# Cells without a cell-type label in the metadata (empty string or NA)
# are kept as their own group, as in the other analyses, and shown as
# "(no label)" (indexing a table by the name "" would return NA).
no_lab <- is.na(ct) | !nzchar(trimws(ct))
cat("Cells without a cell-type label: ", sum(no_lab), "\n", sep = "")
ct[no_lab] <- "(no label)"
keep_ct <- names(which(table(ct) >= MIN_CELLS))
cat("Cell types with >= ", MIN_CELLS, " cells: ", length(keep_ct), " of ",
    length(unique(na.omit(ct))), "\n", sep = "")

mz <- t(sapply(keep_ct, function(x) colMeans(Z[which(ct == x), , drop = FALSE])))
n_ct <- vapply(keep_ct, function(x) sum(ct == x), integer(1))
eps <- sapply(nb$b_CELL, function(b) rev_kw_epsilon2(U[, b], ct)$epsilon2)

ord <- hclust(dist(mz), method = "average")$order
mz <- mz[ord, , drop = FALSE]; n_ct <- n_ct[ord]

out <- data.table(cell_type = rownames(mz), n_cells = n_ct, mz)
fwrite(out, OUT_TSV, sep = "\t")
cat("\nMean standardized loading per cell type:\n")
print(out[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

cols <- hcl.colors(101, "RdBu", rev = TRUE)

draw <- function() {
    m <- pmax(pmin(mz, CLIP), -CLIP)
    nr <- nrow(m); nc <- ncol(m)
    par(mar = c(3, 13, 6.5, 8), xpd = FALSE)
    image(seq_len(nc), seq_len(nr), t(m[nr:1, , drop = FALSE]), col = cols,
          zlim = c(-CLIP, CLIP), axes = FALSE, xlab = "", ylab = "")
    abline(v = seq(0.5, nc + 0.5), h = seq(0.5, nr + 0.5), col = "white", lwd = 0.6)
    box()
    axis(2, at = seq_len(nr), labels = rev(rownames(m)), las = 2, cex.axis = 0.7, tick = FALSE)
    hl <- comps == HIGHLIGHT
    if (any(!hl)) axis(3, at = which(!hl), labels = comps[!hl], las = 2, cex.axis = 0.75,
                       tick = FALSE, line = 2.2, font = 1)
    if (any(hl)) {
        axis(3, at = which(hl), labels = comps[hl], las = 2, cex.axis = 0.75,
             tick = FALSE, line = 2.2, font = 2, col.axis = "firebrick4")
        rect(which(hl) - 0.5, 0.5, which(hl) + 0.5, nr + 0.5, border = "firebrick4", lwd = 2)
    }
    par(xpd = NA)
    text(seq_len(nc), nr + 0.9, labels = sprintf("%d\n%.2f", nb$N, eps), cex = 0.55)
    text(0.3, nr + 0.9, labels = expression(atop("n top-100 cores", "cell-type " * epsilon^2)), adj = c(1, 0.5), cex = 0.55)
    text(nc + 0.6, seq_len(nr), labels = rev(n_ct), adj = c(0, 0.5), cex = 0.6, col = "grey30")
    text(nc + 0.6, nr + 0.9, labels = "n cells", adj = c(0, 0.5), cex = 0.6, col = "grey30")
    par(xpd = FALSE)
    rev_colorbar(cols, c(-CLIP, CLIP), "mean z", at = c(-2, -1, 0, 1, 2),
                 x_from = 1.10, x_to = 1.14)
}

h <- max(5, 0.22 * nrow(mz) + 3)
pdf(OUT_PDF, width = 2.5 + 0.45 * ncol(mz) + 2, height = h, useDingbats = FALSE)
draw(); dev.off()
png(OUT_PNG, width = (2.5 + 0.45 * ncol(mz) + 2) * 300, height = h * 300, res = 300)
draw(); dev.off()

cat("\nSaved:\n", OUT_TSV, "\n", OUT_PDF, "\n", OUT_PNG, "\n", sep = "")
rev_end_log(log_con)
