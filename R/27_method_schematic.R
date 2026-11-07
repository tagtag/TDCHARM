# ============================================================
# 27_method_schematic.R
#
# Revision figure (Reviewer 2, minor comment 1): graphical summary
# of the method, from the four modalities through the Stage-1
# decompositions and the implicit tensor to the Stage-2 HOSVD and
# the downstream analysis. Drawn with base graphics only; no data
# are read. Dimensions in the boxes are those of the manuscript.
#
# Outputs:
#   GSE303006_M23_rev27_method_schematic.pdf
#   GSE303006_M23_rev27_method_schematic.png
# ============================================================

OUT_PDF <- "GSE303006_M23_rev27_method_schematic.pdf"
OUT_PNG <- "GSE303006_M23_rev27_method_schematic.png"

mod_cols <- c(ATAC = "#4C78A8", H3K27ac = "#F58518", RNA = "#54A24B", `3D` = "#B279A2")

box_ <- function(x, y, w, h, label, fill = "white", border = "grey20",
                 cex = 0.8, font = 1, lwd = 1.2) {
    rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill,
         border = border, lwd = lwd)
    text(x, y, label, cex = cex, font = font)
}
arrow_ <- function(x0, y0, x1, y1, col = "grey25") {
    arrows(x0, y0, x1, y1, length = 0.08, lwd = 1.4, col = col)
}
stage_title <- function(x, y, txt) text(x, y, txt, font = 2, cex = 0.95)

draw_schematic <- function() {
    op <- par(mar = c(0.5, 0.5, 0.5, 0.5))
    on.exit(par(op))
    plot(NA, xlim = c(0, 100), ylim = c(0, 62), axes = FALSE, xlab = "", ylab = "",
         asp = NA)

    # --- column 1: four modalities --------------------------------
    stage_title(11, 60, "Input: same 4,258 cells")
    ys <- c(51, 39, 27, 15)
    labs <- c("ATAC\n42,669 enhancers x cells",
              "H3K27ac\n42,669 enhancers x cells",
              "RNA\n16,239 genes x cells",
              "3D (-log distance)\n391,435 bin pairs x cells")
    for (i in 1:4) box_(11, ys[i], 19, 8, labs[i], fill = adjustcolor(mod_cols[i], 0.18),
                        border = mod_cols[i], cex = 0.68)

    # --- column 2: Stage 1 -----------------------------------------
    stage_title(34, 60, "Stage 1: per-modality SVD")
    for (i in 1:4) {
        box_(34, ys[i], 17, 8,
             expression(X[m] %~~% B[m] * V[m]^T),
             fill = "white", border = mod_cols[i], cex = 0.8)
        text(34, ys[i] - 2.9, "K = 20, ||rank-20 part||_F = 1", cex = 0.55, col = "grey30")
        arrow_(20.6, ys[i], 25.4, ys[i], mod_cols[i])
    }

    # --- column 3: mapping to E-P pairs ------------------------------
    stage_title(57, 60, "Map to E-P pairs")
    box_(57, 33, 18, 34,
         paste0("730,969 candidate\nenhancer-promoter pairs\n(enhancer within 1 Mb of TSS)\n\n",
                "each feature copied to its\nd pairs with weight 1/sqrt(d)\n",
                "(d: number of pairs sharing\nthe enhancer, gene or\n20-kb bin pair)\n\n",
                "a_mk(p) = B_m[feature(p),k] / sqrt(d)"),
         fill = "grey97", cex = 0.6)
    for (i in 1:4) arrow_(42.6, ys[i], 47.9, 33 + (2.5 - i) * 4, mod_cols[i])

    # --- column 4: implicit tensor + Stage 2 -------------------------
    stage_title(82, 60, "Stage 2: implicit HOSVD")
    box_(82, 49, 30, 13,
         paste0("Implicit 4-way tensor (never stored)\n",
                "X[p,c,k,m] = a_mk(p) v_mk(c)\n",
                "E-P pair x cell x component x modality\n",
                "730,969 x 4,258 x 20 x 4"),
         fill = "#FFF7E0", border = "#C9A227", cex = 0.62)
    arrow_(66.1, 45, 66.9, 45)
    arrow_(82, 42.3, 82, 38.7)
    box_(82, 31, 30, 15,
         paste0("Factors from small matrices:\n",
                "E-P: eigenvectors of 80 x 80 A^T A\n",
                "cell: SVD of [ ||a_mk|| v_mk ] (4,258 x 80)\n",
                "component: identity (Gram is diagonal)\n",
                "modality: 4 x 4 Gram\n",
                "core G(l1 l2 k m), Tucker ranks 20,20,20,4"),
         fill = "#EEF4FB", border = "#4C78A8", cex = 0.58)
    arrow_(82, 23.3, 82, 19.7)
    box_(82, 12, 30, 15,
         paste0("Downstream analysis\n",
                "1. select cell factor (cell types)\n",
                "2. pick E-P factor through core G\n",
                "3. degree-corrected gene ranking\n",
                "4. GO-BP GSEA\n",
                "5. back-projection to the 4 measurements"),
         fill = "#EEF7EC", border = "#54A24B", cex = 0.58)
}

pdf(OUT_PDF, width = 11, height = 6.5, useDingbats = FALSE)
draw_schematic(); dev.off()
png(OUT_PNG, width = 11 * 300, height = 6.5 * 300, res = 300)
draw_schematic(); dev.off()
cat("Saved:\n", OUT_PDF, "\n", OUT_PNG, "\n", sep = "")
