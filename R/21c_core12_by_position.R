# ============================================================
# 21c_core12_by_position.R
# GSE303006 / GENCODE M23
#
# Follow-up to 21 / 21b (Reviewer 1, Comment 3).
#
# 21b showed that the E-P factor and the cell factor of manuscript
# core #12 are recovered at the SAME positions in every independently
# decomposed two-replicate subset (EP03 |r| ~ 0.97, CELL03 |r| >= 0.91),
# whereas the largest-|core| rule of 21 picked a different, nearly
# equally large core element in the subsets.
#
# Here the core element at the same position as core #12,
#     (EP a, CELL b, K k, MOD r) taken from the full-data top-100 table,
# is analysed in every run. The positions are each run's own ranks, so
# no factor, gene ranking or GO result of the full data is used to
# analyse a subset.
#
# For each run (full data and three subsets):
#   - the core value at that position, its rank among all core
#     elements, and its rank among the elements containing CELL b
#   - |correlation| of the modality-factor column r with the full data
#   - degree-corrected gene score (as in 11), oriented as in 21
#     (cell column: group A higher; modality column: positive sum)
#   - GO-BP GSEA (as in 12; serial)
# and the agreement of the gene rankings and GO results between runs.
#
# Inputs (working directory): outputs of 08-10 and the sub-directories
# rev21_* written by 21.
# Outputs:
#   GSE303006_M23_rev21c_core_by_position.tsv
#   GSE303006_M23_rev21c_gene_scores.tsv.gz
#   GSE303006_M23_rev21c_gene_ranking_agreement.tsv
#   GSE303006_M23_rev21c_GO_agreement.tsv
#   GSE303006_M23_rev21c_GO_top10_per_run.tsv
#   GSE303006_M23_rev21c_GO_all_runs.tsv.gz
#   GSE303006_M23_rev21c_log.txt
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

CORE_RANK <- 12L
SUBSETS <- c("R1R2", "R1R3", "R2R3")
TOP_N <- c(100L, 500L)
FDR <- 0.05
RUN_GSEA <- !identical(Sys.getenv("TDCHARM_RUN_GSEA"), "0") &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE)

CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
EP_FACTOR_FILE <- "GSE303006_M23_stage2_Tucker_EP_factor.rds"
MAPPING_FILE <- "GSE303006_M23_stage2_EP_mapping_degrees.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_CORE <- "GSE303006_M23_rev21c_core_by_position.tsv"
OUT_GENES <- "GSE303006_M23_rev21c_gene_scores.tsv.gz"
OUT_RANK <- "GSE303006_M23_rev21c_gene_ranking_agreement.tsv"
OUT_GO <- "GSE303006_M23_rev21c_GO_agreement.tsv"
OUT_GO_TOP <- "GSE303006_M23_rev21c_GO_top10_per_run.tsv"
OUT_GO_ALL <- "GSE303006_M23_rev21c_GO_all_runs.tsv.gz"
OUT_LOG <- "GSE303006_M23_rev21c_log.txt"

dirs <- c(FULL = ".", setNames(paste0("rev21_", SUBSETS), SUBSETS))
need <- c(TOP_CORE_FILE, META_FILE, MAPPING_FILE,
          as.vector(outer(dirs, c(CORE_FILE, SMALL_FILE, EP_FACTOR_FILE), file.path)))
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "21c_core12_by_position.R")
cat("GSEA: ", if (RUN_GSEA) "run (serial)" else "skipped", "\n", sep = "")

top <- fread(TOP_CORE_FILE)
pos <- top[core_rank == CORE_RANK]
a <- pos$a_EP; b <- pos$b_CELL; k <- pos$k; r <- pos$r_MODAL
cat(sprintf("Position of manuscript core #%d: (EP%02d, CELL%02d, K%02d, MOD%02d)\n\n",
            CORE_RANK, a, b, k, r))

meta <- fread(META_FILE)
groups <- rev_contrast_groups(meta$cell_type)
map <- as.data.table(readRDS(MAPPING_FILE))
U_MODAL_full <- readRDS(file.path(".", SMALL_FILE))$U_MODAL

runs <- list()
for (nm in names(dirs)) {
    d <- dirs[[nm]]
    CORE <- readRDS(file.path(d, CORE_FILE))
    SF <- readRDS(file.path(d, SMALL_FILE))
    UEP <- readRDS(file.path(d, EP_FACTOR_FILE))
    mp <- as.data.table(readRDS(file.path(d, MAPPING_FILE)))
    stopifnot(identical(mp$gene_name, map$gene_name))
    cells <- rownames(SF$U_CELL)
    ct <- meta$cell_type[match(cells, meta$cell)]
    inA <- ct %in% groups$A; inB <- ct %in% groups$B

    G <- CORE[a, b, k, r]
    s_b <- if (mean(SF$U_CELL[inA, b]) >= mean(SF$U_CELL[inB, b])) 1 else -1
    s_r <- if (sum(SF$U_MODAL[, r]) >= 0) 1 else -1
    G_or <- G * s_b * s_r
    gene <- rev_gene_scores(mp, G_or * UEP[, a])
    gene[, run := nm]

    slice_b <- abs(as.numeric(CORE[, b, , ]))
    runs[[nm]] <- list(
        gene = gene,
        row = data.table(
            run = nm, n_cells = length(cells),
            EP = a, CELL = b, K = k, MOD = r,
            core_value_raw = G, core_value_oriented = G_or,
            rank_among_all_core = sum(abs(as.numeric(CORE)) > abs(G)) + 1L,
            rank_among_CELL_elements = sum(slice_b > abs(G)) + 1L,
            ratio_to_largest_CELL_element = abs(G) / max(slice_b),
            modality_column_abs_r_with_full = abs(cor(SF$U_MODAL[, r], U_MODAL_full[, r]))))
    rm(UEP); invisible(gc())
}

core_dt <- rbindlist(lapply(runs, `[[`, "row"))
cat("Core element at the position of core #12 in every run:\n")
print(core_dt[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)])

pairs <- combn(names(runs), 2, simplify = FALSE)
rank_agree <- rbindlist(lapply(pairs, function(p)
    rev_compare_rankings(runs[[p[1]]]$gene, runs[[p[2]]]$gene, p[1], p[2], TOP_N)))
cat("\nGene-ranking agreement:\n")
print(rank_agree[, .(run_x, run_y, spearman = round(spearman, 4), side, top_n,
                     overlap, expected = round(expected, 2),
                     hypergeom_p = signif(hypergeom_p, 3))])

go_agree <- data.table(); go_top <- data.table(); go_all <- data.table()
if (RUN_GSEA) {
    for (nm in names(runs)) {
        cat("  GSEA for ", nm, "...\n", sep = "")
        g <- rev_run_gsea(runs[[nm]]$gene)
        if (nrow(g)) g[, run := nm]
        runs[[nm]]$go <- g
    }
    go_agree <- rbindlist(lapply(pairs, function(p)
        rev_compare_go(runs[[p[1]]]$go, runs[[p[2]]]$go, p[1], p[2], FDR)))
    go_top <- rbindlist(lapply(names(runs), function(nm) rev_go_top_terms(runs[[nm]]$go, nm, 10L)))
    go_all <- rbindlist(lapply(runs, `[[`, "go"), fill = TRUE)
    cat("\nGO agreement (FDR <= ", FDR, "):\n", sep = "")
    print(go_agree)
    cat("\nTop 10 GO terms per side and run (no omission):\n")
    print(go_top[, .(run, side, rank, ID, Description = substr(Description, 1, 60),
                     NES = round(NES, 2), p.adjust = signif(p.adjust, 3))])
}

fwrite(core_dt, OUT_CORE, sep = "\t")
fwrite(rbindlist(lapply(runs, `[[`, "gene")), OUT_GENES, sep = "\t")
fwrite(rank_agree, OUT_RANK, sep = "\t")
if (RUN_GSEA) {
    fwrite(go_agree, OUT_GO, sep = "\t")
    fwrite(go_top, OUT_GO_TOP, sep = "\t")
    fwrite(go_all[, lapply(.SD, function(x) if (is.list(x)) vapply(x, paste, "", collapse = "/") else x)],
           OUT_GO_ALL, sep = "\t")
}
cat("\nSaved outputs with prefix GSE303006_M23_rev21c_.\n")
rev_end_log(log_con)
