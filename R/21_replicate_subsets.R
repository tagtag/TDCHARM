# ============================================================
# 21_replicate_subsets.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 1, Comment 3).
#
# Is the reported program driven by any single replicate?
#
# The whole decomposition is repeated INDEPENDENTLY on each of the
# three two-replicate subsets (R1+R2, R1+R3, R2+R3), using exactly
# the released Stage-1 / Stage-2 scripts 05-08 (sourced unchanged,
# each in its own sub-directory, with the cell list restricted to
# the subset). No factor, core, gene ranking or GO result from the
# full analysis is used to analyse a subset.
#
# In every run (full data and the three subsets) the same
# PRE-SPECIFIED rule selects the core element (revision_lib.R,
# rev_select_core): the cell component whose 50 most positive and
# 50 most negative cells are most enriched for the two cell types
# characterizing core #12 (Inh Ndnf/Lamp5 versus Ex L3/4 IT), then
# the largest |core| element containing that cell component. The
# rule depends only on the cell-type annotation, not on the full
# analysis.
#
# For each run:
#   - degree-corrected gene score (as in 11)
#   - preranked GO-BP GSEA (as in 12; skipped if clusterProfiler is
#     unavailable or TDCHARM_RUN_GSEA=0)
#   - cell-type and replicate epsilon^2 of the selected cell
#     component, and of all cell components in the top-100 core
#     elements (as in 10/17)
# Final results are then compared between every pair of runs:
#   - Spearman correlation of gene scores, top-N gene overlap
#     (hypergeometric P)
#   - overlap / Jaccard of FDR <= 0.05 GO terms per side, Spearman
#     correlation of NES
#   - correlation of the selected cell loadings on shared cells
#
# Subset runs are cached: if a sub-directory already contains the
# Stage-2 core it is reused unless TDCHARM_RECOMPUTE=1.
#
# Inputs (working directory): preprocessing outputs (03/04), the
# full-data outputs of 06-11, and the metadata table of 10.
# Outputs:
#   rev21_<subset>/            Stage-1/Stage-2 outputs per subset
#   GSE303006_M23_rev21_selection.tsv
#   GSE303006_M23_rev21_selection_component_scores.tsv
#   GSE303006_M23_rev21_selection_core_candidates.tsv
#   GSE303006_M23_rev21_gene_scores.tsv.gz
#   GSE303006_M23_rev21_gene_ranking_agreement.tsv
#   GSE303006_M23_rev21_GO_agreement.tsv
#   GSE303006_M23_rev21_GO_top10_per_run.tsv
#   GSE303006_M23_rev21_GO_all_runs.tsv.gz
#   GSE303006_M23_rev21_cell_side.tsv
#   GSE303006_M23_rev21_top100_cell_components.tsv
#   GSE303006_M23_rev21_log.txt
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

REPO <- rev_find_repo()
WORKDIR <- normalizePath(getwd())

SUBSETS <- list(
    R1R2 = c("R1", "R2"),
    R1R3 = c("R1", "R3"),
    R2R3 = c("R2", "R3")
)

TOP_N <- c(100L, 500L)
TOP_CORE <- 100L
FDR <- 0.05
MANUSCRIPT_CORE_RANK <- 12L
RECOMPUTE <- identical(Sys.getenv("TDCHARM_RECOMPUTE"), "1")
RUN_GSEA <- !identical(Sys.getenv("TDCHARM_RUN_GSEA"), "0") &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE)

# Files read by 05 (linked into each sub-directory)
STAGE1_INPUTS <- c(
    "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds",
    "GSE303006_M23_brain_final_cell_QC.tsv",
    "GSE303006_M23_brain_RNA_logCP10K_gene_by_cell.rds",
    "GSE303006_M23_brain_ATAC_enhancer_by_cell_AH3D.rds",
    "GSE303006_M23_brain_H3K27ac_enhancer_by_cell_AH3D.rds",
    "GSE303006_M23_brain_3D_binpair_by_cell.h5"
)
# 05 reads its cell list from this file name; in a sub-directory it
# contains only the cells of the subset.
CELL_FILE <- "GSE303006_M23_brain_final_4258_cells.rds"

META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"
MAPPING_FILE <- "GSE303006_M23_stage2_EP_mapping_degrees.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
EP_FACTOR_FILE <- "GSE303006_M23_stage2_Tucker_EP_factor.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"
GENE12_FILE <- "GSE303006_M23_core12_all16239_gene_rank_degree_corrected.tsv"
GO12_FILE <- "GSE303006_M23_core12_GO_BP_GSEA_all.tsv"

OUT_SEL <- "GSE303006_M23_rev21_selection.tsv"
OUT_SEL_COMP <- "GSE303006_M23_rev21_selection_component_scores.tsv"
OUT_CAND <- "GSE303006_M23_rev21_selection_core_candidates.tsv"
OUT_GENES <- "GSE303006_M23_rev21_gene_scores.tsv.gz"
OUT_RANK_AGREE <- "GSE303006_M23_rev21_gene_ranking_agreement.tsv"
OUT_GO_AGREE <- "GSE303006_M23_rev21_GO_agreement.tsv"
OUT_GO_TOP <- "GSE303006_M23_rev21_GO_top10_per_run.tsv"
OUT_GO_ALL <- "GSE303006_M23_rev21_GO_all_runs.tsv.gz"
OUT_CELL <- "GSE303006_M23_rev21_cell_side.tsv"
OUT_TOP100 <- "GSE303006_M23_rev21_top100_cell_components.tsv"
OUT_LOG <- "GSE303006_M23_rev21_log.txt"

need <- c(STAGE1_INPUTS, CELL_FILE, META_FILE, MAPPING_FILE, CORE_FILE,
          SMALL_FILE, EP_FACTOR_FILE, TOP_CORE_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "21_replicate_subsets.R")
cat("Repository: ", REPO, "\n", sep = "")
cat("GSEA: ", if (RUN_GSEA) "run" else "skipped", "\n", sep = "")
cat("Recompute cached subsets: ", RECOMPUTE, "\n\n", sep = "")

meta <- fread(META_FILE)
stopifnot(all(c("cell", "cell_type", "replicate") %in% names(meta)))
groups <- rev_contrast_groups(meta$cell_type)
cat("Pre-specified contrast:\n  A = ", paste(groups$A, collapse = "; "),
    "\n  B = ", paste(groups$B, collapse = "; "), "\n\n", sep = "")

full_cells <- as.character(readRDS(CELL_FILE))
map <- as.data.table(readRDS(MAPPING_FILE))

# ------------------------------------------------------------
# 1. Independent Stage-1/Stage-2 runs on each subset
# ------------------------------------------------------------

run_pipeline_in <- function(subdir, cells) {
    dir.create(subdir, showWarnings = FALSE)
    for (f in STAGE1_INPUTS) {
        dst <- file.path(subdir, f)
        if (!file.exists(dst)) file.symlink(file.path(WORKDIR, f), dst)
    }
    saveRDS(cells, file.path(subdir, CELL_FILE))
    old <- setwd(subdir)
    on.exit(setwd(old), add = TRUE)
    for (s in c("05_stage1_independent_K20.R", "06_sign_only_stage1_K20.R",
                "07_prepare_stage2_implicit.R", "08_stage2_implicit_HOSVD.R")) {
        cat("  [", basename(subdir), "] running ", s, "\n", sep = "")
        out <- capture.output(
            source(file.path(REPO, "R", s), local = new.env(parent = globalenv())),
            type = "output")
        writeLines(out, paste0("log_", sub("\\.R$", "", s), ".txt"))
    }
    invisible(NULL)
}

for (nm in names(SUBSETS)) {
    subdir <- file.path(WORKDIR, paste0("rev21_", nm))
    cells <- full_cells[substr(full_cells, 1, 2) %in% SUBSETS[[nm]]]
    cat("Subset ", nm, ": ", length(cells), " cells\n", sep = "")
    if (!RECOMPUTE && file.exists(file.path(subdir, CORE_FILE))) {
        cat("  cached Stage-2 output found; reusing\n")
        next
    }
    run_pipeline_in(subdir, cells)
}
cat("\n")

# ------------------------------------------------------------
# 2. Analyse each run with the pre-specified rule
# ------------------------------------------------------------

analyse_run <- function(label, dir) {
    CORE <- readRDS(file.path(dir, CORE_FILE))
    SF <- readRDS(file.path(dir, SMALL_FILE))
    UEP <- readRDS(file.path(dir, EP_FACTOR_FILE))
    map_r <- as.data.table(readRDS(file.path(dir, MAPPING_FILE)))
    stopifnot(identical(map_r$enhancer_id, map$enhancer_id),
              identical(map_r$gene_name, map$gene_name))
    cells <- rownames(SF$U_CELL)
    mm <- meta[match(cells, cell)]
    sel <- rev_select_core(SF$U_CELL, CORE, SF$U_MODAL, mm$cell_type, groups)
    gene <- rev_gene_scores(map_r, sel$G_oriented * UEP[, sel$a])
    gene[, run := label]

    x <- sel$cell_loading_oriented
    kw_ct <- rev_kw_epsilon2(x, mm$cell_type)
    kw_rep <- rev_kw_epsilon2(x, mm$replicate)

    # top-100 core cell components (as in 10 / 17)
    ord <- order(abs(as.numeric(CORE)), decreasing = TRUE)[seq_len(TOP_CORE)]
    bs <- arrayInd(ord, dim(CORE))[, 2]
    tb <- table(bs)
    top100 <- rbindlist(lapply(as.integer(names(tb)), function(b) {
        e1 <- rev_kw_epsilon2(SF$U_CELL[, b], mm$cell_type)
        e2 <- rev_kw_epsilon2(SF$U_CELL[, b], mm$replicate)
        data.table(run = label, CELL_component = sprintf("CELL%02d", b),
                   n_top100_cores = as.integer(tb[as.character(b)]),
                   celltype_epsilon2 = e1$epsilon2, replicate_epsilon2 = e2$epsilon2)
    }))

    go <- data.table()
    if (RUN_GSEA) {
        cat("  GSEA for ", label, "...\n", sep = "")
        go <- rev_run_gsea(gene)
        if (nrow(go)) go[, run := label]
    }

    list(label = label, sel = sel, gene = gene, go = go,
         loading = setNames(x, cells), top100 = top100,
         sel_row = data.table(
             run = label, n_cells = length(cells),
             EP = sel$a, CELL = sel$b, K = sel$k, MOD = sel$r,
             core_rank_of_selected = sel$core_rank,
             core_value_raw = sel$G_raw, core_value_oriented = sel$G_oriented,
             selection_score = sel$component_scores$score[1],
             runner_up_score = sel$component_scores$score[2],
             abs_rank_biserial_A_vs_B = abs(sel$component_scores$rank_biserial_A_vs_B[1]),
             celltype_epsilon2 = kw_ct$epsilon2, celltype_KW_p = kw_ct$p,
             replicate_epsilon2 = kw_rep$epsilon2, replicate_KW_p = kw_rep$p),
         comp_scores = cbind(run = label, sel$component_scores),
         candidates = cbind(run = label, sel$core_candidates))
}

runs <- list()
runs$FULL <- analyse_run("FULL", WORKDIR)
for (nm in names(SUBSETS)) {
    runs[[nm]] <- analyse_run(nm, file.path(WORKDIR, paste0("rev21_", nm)))
}

sel_dt <- rbindlist(lapply(runs, `[[`, "sel_row"))
cat("Selected core elements (pre-specified rule):\n")
print(sel_dt[, .(run, n_cells, EP, CELL, K, MOD, core_rank_of_selected,
                 selection_score = round(selection_score, 3),
                 runner_up_score = round(runner_up_score, 3),
                 celltype_eps2 = round(celltype_epsilon2, 3),
                 replicate_eps2 = round(replicate_epsilon2, 3))])

cat("\nLargest core elements containing the selected cell component\n",
    "(candidate 1 is selected; a ratio close to 1 for candidate 2 means a near tie):\n", sep = "")
cand_dt <- rbindlist(lapply(runs, `[[`, "candidates"))
print(cand_dt[, .(run, candidate, EP, CELL, K, MOD, abs_core = signif(abs_core, 4),
                  ratio_to_selected = round(ratio_to_selected, 3))])

# Relation of the rule-selected full-data core to manuscript core #12
top <- fread(TOP_CORE_FILE)
c12 <- top[core_rank == MANUSCRIPT_CORE_RANK]
f <- runs$FULL$sel
same12 <- f$a == c12$a_EP && f$b == c12$b_CELL && f$k == c12$k && f$r == c12$r_MODAL
cat(sprintf("\nManuscript core #%d = (EP%02d, CELL%02d, K%02d, MOD%02d); rule-selected full-data core = (EP%02d, CELL%02d, K%02d, MOD%02d): %s\n",
            MANUSCRIPT_CORE_RANK, c12$a_EP, c12$b_CELL, c12$k, c12$r_MODAL,
            f$a, f$b, f$k, f$r, if (same12) "IDENTICAL" else "DIFFERENT"))
if (file.exists(GENE12_FILE)) {
    g12 <- fread(GENE12_FILE)[, .(gene_name, s12 = gene_score_degree_corrected)]
    mm <- merge(runs$FULL$gene, g12, by = "gene_name")
    cat(sprintf("Spearman(rule-selected full-data gene score, manuscript core #%d gene score) = %.4f\n",
                MANUSCRIPT_CORE_RANK, cor(mm$gene_score, mm$s12, method = "spearman")))
}
if (!same12) warning("Rule-selected full-data core differs from the manuscript core; see log.")

# ------------------------------------------------------------
# 3. Agreement between runs
# ------------------------------------------------------------

pairs <- combn(names(runs), 2, simplify = FALSE)

rank_agree <- rbindlist(lapply(pairs, function(p)
    rev_compare_rankings(runs[[p[1]]]$gene, runs[[p[2]]]$gene, p[1], p[2], TOP_N)))
cat("\nGene-ranking agreement between runs:\n")
print(rank_agree[, .(run_x, run_y, spearman = round(spearman, 4), side, top_n,
                     overlap, expected = round(expected, 2),
                     hypergeom_p = signif(hypergeom_p, 3))])

go_agree <- data.table()
go_top <- data.table()
go_all <- data.table()
if (RUN_GSEA) {
    go_agree <- rbindlist(lapply(pairs, function(p)
        rev_compare_go(runs[[p[1]]]$go, runs[[p[2]]]$go, p[1], p[2], FDR)))
    go_top <- rbindlist(lapply(runs, function(r) rev_go_top_terms(r$go, r$label, 10L)))
    go_all <- rbindlist(lapply(runs, `[[`, "go"), fill = TRUE)
    cat("\nGO agreement between runs (FDR <= ", FDR, "):\n", sep = "")
    print(go_agree)
    cat("\nTop 10 GO terms per side and run (no omission):\n")
    print(go_top[, .(run, side, rank, ID, Description = substr(Description, 1, 60),
                     NES = round(NES, 2), p.adjust = signif(p.adjust, 3))])
}

cell_side <- rbindlist(lapply(pairs, function(p) {
    lx <- runs[[p[1]]]$loading; ly <- runs[[p[2]]]$loading
    common <- intersect(names(lx), names(ly))
    data.table(run_x = p[1], run_y = p[2], n_shared_cells = length(common),
               pearson_selected_cell_loading = cor(lx[common], ly[common]),
               spearman_selected_cell_loading = cor(lx[common], ly[common], method = "spearman"))
}))
cat("\nSelected cell component: agreement on shared cells:\n")
print(cell_side)

top100 <- rbindlist(lapply(runs, `[[`, "top100"))
cat("\nTop-100 core cell components per run: cell type vs replicate epsilon^2\n")
print(top100[, .(n_components = .N,
                 n_celltype_gt_replicate = sum(celltype_epsilon2 > replicate_epsilon2),
                 median_celltype_eps2 = median(celltype_epsilon2),
                 median_replicate_eps2 = median(replicate_epsilon2)), by = run])

# ------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------

fwrite(sel_dt, OUT_SEL, sep = "\t")
fwrite(rbindlist(lapply(runs, `[[`, "comp_scores")), OUT_SEL_COMP, sep = "\t")
fwrite(cand_dt, OUT_CAND, sep = "\t")
fwrite(rbindlist(lapply(runs, `[[`, "gene")), OUT_GENES, sep = "\t")
fwrite(rank_agree, OUT_RANK_AGREE, sep = "\t")
fwrite(cell_side, OUT_CELL, sep = "\t")
fwrite(top100, OUT_TOP100, sep = "\t")
if (RUN_GSEA) {
    fwrite(go_agree, OUT_GO_AGREE, sep = "\t")
    fwrite(go_top, OUT_GO_TOP, sep = "\t")
    fwrite(go_all[, lapply(.SD, function(x) if (is.list(x)) vapply(x, paste, "", collapse = "/") else x)],
           OUT_GO_ALL, sep = "\t")
}

cat("\nSaved outputs with prefix GSE303006_M23_rev21_ and sub-directories rev21_*.\n")
rev_end_log(log_con)
