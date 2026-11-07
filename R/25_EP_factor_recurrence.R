# ============================================================
# 25_EP_factor_recurrence.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 2, Comment 4: pair-mode factors).
#
# 1. How often does each E-P factor occur among the 100 largest
#    core elements (count and summed |core|)?
# 2. Which (cell factor, E-P factor) combinations recur among them?
#    (20 x 20 table of counts and summed |core|, with heatmaps.)
#    Each cell factor is described by the cell types enriched
#    among its 50 most positive and 50 most negative cells.
# 3. The N_EP most frequent E-P factors are characterized through
#    their representative core element (the largest |core| element
#    with that E-P factor):
#      - degree-corrected gene ranking (as in 11), top genes
#      - preranked GO-BP GSEA (as in 12; skipped if clusterProfiler
#        is unavailable or TDCHARM_RUN_GSEA=0), top 10 terms per
#        side without omission
#      - enhancer-TSS distance of the 1,000 highest-|score| pairs
#        compared with all candidate pairs
#
# Inputs (working directory): outputs of 03, 08, 09, 10
# Outputs:
#   GSE303006_M23_rev25_EP_factor_frequency.tsv
#   GSE303006_M23_rev25_EP_by_CELL_counts.tsv
#   GSE303006_M23_rev25_EP_by_CELL_abs_core.tsv
#   GSE303006_M23_rev25_EP_by_CELL_heatmap.pdf
#   GSE303006_M23_rev25_cell_component_ends.tsv
#   GSE303006_M23_rev25_recurrent_EP_top_genes.tsv
#   GSE303006_M23_rev25_recurrent_EP_GO_top10.tsv
#   GSE303006_M23_rev25_recurrent_EP_GO_all.tsv.gz
#   GSE303006_M23_rev25_recurrent_EP_distance.tsv
#   GSE303006_M23_rev25_log.txt
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

# ------------------------------------------------------------
# 0. Settings / files
# ------------------------------------------------------------

TOP_CORE <- 100L
N_EP <- as.integer(Sys.getenv("TDCHARM_N_RECURRENT_EP", unset = "3"))
TOP_CELLS <- 50L
TOP_GENES <- 20L
TOP_PAIRS_DIST <- 1000L
RUN_GSEA <- !identical(Sys.getenv("TDCHARM_RUN_GSEA"), "0") &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE)

EP_FILE <- "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"
MAPPING_FILE <- "GSE303006_M23_stage2_EP_mapping_degrees.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
EP_FACTOR_FILE <- "GSE303006_M23_stage2_Tucker_EP_factor.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_FREQ <- "GSE303006_M23_rev25_EP_factor_frequency.tsv"
OUT_CNT <- "GSE303006_M23_rev25_EP_by_CELL_counts.tsv"
OUT_ABS <- "GSE303006_M23_rev25_EP_by_CELL_abs_core.tsv"
OUT_PDF <- "GSE303006_M23_rev25_EP_by_CELL_heatmap.pdf"
OUT_ENDS <- "GSE303006_M23_rev25_cell_component_ends.tsv"
OUT_GENES <- "GSE303006_M23_rev25_recurrent_EP_top_genes.tsv"
OUT_GO_TOP <- "GSE303006_M23_rev25_recurrent_EP_GO_top10.tsv"
OUT_GO_ALL <- "GSE303006_M23_rev25_recurrent_EP_GO_all.tsv.gz"
OUT_DIST <- "GSE303006_M23_rev25_recurrent_EP_distance.tsv"
OUT_LOG <- "GSE303006_M23_rev25_log.txt"

need <- c(EP_FILE, MAPPING_FILE, CORE_FILE, SMALL_FILE, EP_FACTOR_FILE,
          TOP_CORE_FILE, META_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "25_EP_factor_recurrence.R")
cat("GSEA: ", if (RUN_GSEA) "run" else "skipped", "\n\n", sep = "")

top <- fread(TOP_CORE_FILE)[core_rank <= TOP_CORE]
CORE <- readRDS(CORE_FILE)
SF <- readRDS(SMALL_FILE)
map <- as.data.table(readRDS(MAPPING_FILE))
ep <- as.data.table(readRDS(EP_FILE))
stopifnot(identical(ep$enhancer_id, map$enhancer_id), identical(ep$gene_name, map$gene_name))
R_EP <- dim(CORE)[1]; R_CELL <- dim(CORE)[2]

# ------------------------------------------------------------
# 1. E-P factor frequency
# ------------------------------------------------------------

freq <- top[, .(n_top_cores = .N, sum_abs_core = sum(abs_core),
                best_core_rank = min(core_rank),
                cell_partners = paste(sort(unique(sprintf("CELL%02d", b_CELL))), collapse = ",")),
            by = .(a_EP)]
freq[, EP_component := sprintf("EP%02d", a_EP)]
setorder(freq, -n_top_cores, -sum_abs_core)
cat("E-P factors among the top-", TOP_CORE, " core elements:\n", sep = "")
print(freq[, .(EP_component, n_top_cores, sum_abs_core = round(sum_abs_core, 4),
               best_core_rank, cell_partners)])

# ------------------------------------------------------------
# 2. (cell, E-P) co-occurrence
# ------------------------------------------------------------

cnt <- matrix(0L, R_EP, R_CELL, dimnames = list(sprintf("EP%02d", seq_len(R_EP)),
                                                sprintf("CELL%02d", seq_len(R_CELL))))
absm <- matrix(0, R_EP, R_CELL, dimnames = dimnames(cnt))
for (i in seq_len(nrow(top))) {
    cnt[top$a_EP[i], top$b_CELL[i]] <- cnt[top$a_EP[i], top$b_CELL[i]] + 1L
    absm[top$a_EP[i], top$b_CELL[i]] <- absm[top$a_EP[i], top$b_CELL[i]] + top$abs_core[i]
}
pairs <- top[, .(n_top_cores = .N, sum_abs_core = sum(abs_core),
                 best_core_rank = min(core_rank)), by = .(a_EP, b_CELL)]
setorder(pairs, -n_top_cores, -sum_abs_core)
cat("\nRecurring (E-P factor, cell factor) combinations (>= 2 of the top-", TOP_CORE, "):\n", sep = "")
print(pairs[n_top_cores >= 2L, .(EP = sprintf("EP%02d", a_EP), CELL = sprintf("CELL%02d", b_CELL),
                                 n_top_cores, sum_abs_core = round(sum_abs_core, 4), best_core_rank)])

# Cell factor ends: cell types among the 50 most positive / negative cells
meta <- fread(META_FILE)[match(rownames(SF$U_CELL), cell)]
ends <- rbindlist(lapply(seq_len(R_CELL), function(b) {
    x <- SF$U_CELL[, b]
    o <- order(x, decreasing = TRUE)
    desc <- function(ii) {
        tt <- sort(table(meta$cell_type[ii]), decreasing = TRUE)
        paste(sprintf("%s (%d)", names(tt)[seq_len(min(3, length(tt)))],
                      as.integer(tt)[seq_len(min(3, length(tt)))]), collapse = "; ")
    }
    data.table(CELL_component = sprintf("CELL%02d", b),
               positive_end_top3 = desc(o[seq_len(TOP_CELLS)]),
               negative_end_top3 = desc(rev(o)[seq_len(TOP_CELLS)]),
               celltype_epsilon2 = rev_kw_epsilon2(x, meta$cell_type)$epsilon2)
}))
cat("\nCell factors in the top-", TOP_CORE, " core elements (cell types at the two ends):\n", sep = "")
print(ends[CELL_component %in% sprintf("CELL%02d", unique(top$b_CELL))])

# heatmaps
draw_heat <- function(mat, main, legend_title) {
    used_r <- which(rowSums(mat) > 0); used_c <- which(colSums(mat) > 0)
    m <- mat[used_r, used_c, drop = FALSE]
    cols <- hcl.colors(100, "YlOrRd", rev = TRUE)
    zlim <- c(0, max(m))
    par(mar = c(5, 6, 4, 7))
    image(seq_len(ncol(m)), seq_len(nrow(m)), t(m[nrow(m):1, , drop = FALSE]),
          col = cols, zlim = zlim, axes = FALSE, xlab = "Cell factor", ylab = "", main = main)
    axis(1, seq_len(ncol(m)), colnames(m), las = 2, cex.axis = 0.8)
    axis(2, seq_len(nrow(m)), rev(rownames(m)), las = 2, cex.axis = 0.8)
    mtext("E-P factor", side = 2, line = 4.5)
    for (i in seq_len(nrow(m))) for (j in seq_len(ncol(m))) {
        v <- m[nrow(m) - i + 1, j]
        if (v > 0) text(j, i, if (is.integer(v)) v else sprintf("%.2f", v), cex = 0.6)
    }
    box()
    rev_colorbar(cols, zlim, legend_title)
}
pdf(OUT_PDF, width = 8, height = 7)
draw_heat(cnt, sprintf("Top-%d core elements: (E-P, cell) factor counts", TOP_CORE), "count")
draw_heat(absm, sprintf("Top-%d core elements: summed |core|", TOP_CORE), "sum |core|")
dev.off()

# ------------------------------------------------------------
# 3. Characterize the most frequent E-P factors
# ------------------------------------------------------------

UEP <- readRDS(EP_FACTOR_FILE)
rec <- head(freq, N_EP)
gene_rows <- list(); go_rows <- list(); go_all <- list(); dist_rows <- list()
for (i in seq_len(nrow(rec))) {
    a <- rec$a_EP[i]
    z <- top[a_EP == a][which.min(core_rank)]
    G <- CORE[z$a_EP, z$b_CELL, z$k, z$r_MODAL]
    s <- G * UEP[, a]
    lab <- sprintf("EP%02d (core rank %d: CELL%02d, K%02d, MOD%02d)", a, z$core_rank,
                   z$b_CELL, z$k, z$r_MODAL)
    cat("\n", lab, "\n", sep = "")
    gene <- rev_gene_scores(map, s)
    gtop <- rbind(gene[order(-gene_score)][seq_len(TOP_GENES)][, side := "positive"],
                  gene[order(gene_score)][seq_len(TOP_GENES)][, side := "negative"])
    gtop[, `:=`(EP_component = sprintf("EP%02d", a), representative_core_rank = z$core_rank,
                rank = seq_len(.N)), by = side]
    gene_rows[[i]] <- gtop
    cat("  top positive genes: ", paste(gtop[side == "positive", gene_name], collapse = ", "), "\n", sep = "")
    cat("  top negative genes: ", paste(gtop[side == "negative", gene_name], collapse = ", "), "\n", sep = "")

    o <- order(abs(s), decreasing = TRUE)[seq_len(min(TOP_PAIRS_DIST, length(s)))]
    wt <- suppressWarnings(wilcox.test(ep$abs_distance[o], ep$abs_distance))
    dist_rows[[i]] <- data.table(
        EP_component = sprintf("EP%02d", a), representative_core_rank = z$core_rank,
        n_top_pairs = length(o),
        median_distance_top = median(ep$abs_distance[o]),
        median_distance_all = median(ep$abs_distance),
        wilcoxon_p = wt$p.value)

    if (RUN_GSEA) {
        go <- rev_run_gsea(gene)
        lab2 <- sprintf("EP%02d", a)
        go_rows[[i]] <- rev_go_top_terms(go, lab2, 10L)
        if (nrow(go)) go_all[[i]] <- go[, run := lab2]
        cat("  GO terms at FDR <= 0.05: positive ", go[p.adjust <= 0.05 & NES > 0, .N],
            ", negative ", go[p.adjust <= 0.05 & NES < 0, .N], "\n", sep = "")
    }
}
dist_dt <- rbindlist(dist_rows)
cat("\nEnhancer-TSS distance of the top-|score| pairs:\n")
print(dist_dt)
if (RUN_GSEA && length(go_rows)) {
    cat("\nTop 10 GO terms per side (no omission):\n")
    print(rbindlist(go_rows)[, .(run, side, rank, ID, Description = substr(Description, 1, 60),
                                 NES = round(NES, 2), p.adjust = signif(p.adjust, 3))])
}

# ------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------

fwrite(freq, OUT_FREQ, sep = "\t")
fwrite(as.data.table(cnt, keep.rownames = "EP_component"), OUT_CNT, sep = "\t")
fwrite(as.data.table(absm, keep.rownames = "EP_component"), OUT_ABS, sep = "\t")
fwrite(ends, OUT_ENDS, sep = "\t")
fwrite(rbindlist(gene_rows), OUT_GENES, sep = "\t")
fwrite(dist_dt, OUT_DIST, sep = "\t")
if (RUN_GSEA) {
    fwrite(rbindlist(go_rows), OUT_GO_TOP, sep = "\t")
    ga <- rbindlist(go_all, fill = TRUE)
    if (nrow(ga)) fwrite(ga[, lapply(.SD, function(x) if (is.list(x)) vapply(x, paste, "", collapse = "/") else x)],
                         OUT_GO_ALL, sep = "\t")
}
cat("\nSaved outputs with prefix GSE303006_M23_rev25_.\n")
rev_end_log(log_con)
