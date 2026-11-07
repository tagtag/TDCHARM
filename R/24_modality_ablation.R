# ============================================================
# 24_modality_ablation.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 1, Comment 1: individual-modality and
# modality-ablation variants).
#
# Stage 1 is kept as it is (each modality is decomposed
# independently, so removing a modality does not change the others).
# Stage 2 is recomputed with the implicit HOSVD (revision_lib.R; the
# same algebra as 08, checked against 08 for the four-modality case)
# for
#   - every single modality            (4 variants)
#   - every leave-one-modality-out set (4 variants)
#   - all four modalities              (reference)
#
# In every variant the same pre-specified rule as in 21 selects the
# core element (cell component best separating Inh Ndnf/Lamp5 from
# Ex L3/4 IT at the extremes; then the largest |core| element with
# that cell component). Reported per variant:
#   - selected core, its cell-type and replicate epsilon^2
#   - cell-type / replicate epsilon^2 of the cell components in the
#     top-100 core elements
#   - agreement of the degree-corrected gene ranking with the
#     four-modality result (Spearman, top-N overlap)
#   - GO-BP GSEA and its agreement with the four-modality result
#     (skipped if clusterProfiler is unavailable or
#     TDCHARM_RUN_GSEA=0)
#
# Inputs (working directory): outputs of 06-10
# Outputs:
#   GSE303006_M23_rev24_variants.tsv
#   GSE303006_M23_rev24_gene_ranking_agreement.tsv
#   GSE303006_M23_rev24_GO_agreement.tsv
#   GSE303006_M23_rev24_GO_top10_per_variant.tsv
#   GSE303006_M23_rev24_top100_cell_components.tsv
#   GSE303006_M23_rev24_gene_scores.tsv.gz
#   GSE303006_M23_rev24_log.txt
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

TOP_N <- c(100L, 500L)
TOP_CORE <- 100L
FDR <- 0.05
RUN_GSEA <- !identical(Sys.getenv("TDCHARM_RUN_GSEA"), "0") &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE)

MANIFEST_FILE <- "GSE303006_M23_stage2_implicit_manifest.rds"
MAPPING_FILE <- "GSE303006_M23_stage2_EP_mapping_degrees.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_VAR <- "GSE303006_M23_rev24_variants.tsv"
OUT_RANK <- "GSE303006_M23_rev24_gene_ranking_agreement.tsv"
OUT_GO <- "GSE303006_M23_rev24_GO_agreement.tsv"
OUT_GO_TOP <- "GSE303006_M23_rev24_GO_top10_per_variant.tsv"
OUT_TOP100 <- "GSE303006_M23_rev24_top100_cell_components.tsv"
OUT_GENES <- "GSE303006_M23_rev24_gene_scores.tsv.gz"
OUT_LOG <- "GSE303006_M23_rev24_log.txt"

need <- c(MANIFEST_FILE, MAPPING_FILE, CORE_FILE, META_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "24_modality_ablation.R")
cat("GSEA: ", if (RUN_GSEA) "run" else "skipped", "\n\n", sep = "")

mods <- c("ATAC", "H3K27ac", "RNA", "3D")
M <- readRDS(MANIFEST_FILE)
map <- as.data.table(readRDS(MAPPING_FILE))
S <- lapply(M$stage1_files, readRDS)[mods]
cells <- rownames(S[[1]]$cell_loadings)
meta <- fread(META_FILE)[match(cells, cell)]
groups <- rev_contrast_groups(meta$cell_type)

variants <- c(list(ALL4 = mods),
              setNames(lapply(mods, function(m) setdiff(mods, m)), paste0("minus_", mods)),
              setNames(as.list(mods), paste0("only_", mods)))

# ------------------------------------------------------------
# 1. Fit every variant
# ------------------------------------------------------------

res <- list()
for (vn in names(variants)) {
    vm <- variants[[vn]]
    fit <- rev_stage2_hosvd(S[vm], map)
    if (vn == "ALL4") {
        C08 <- readRDS(CORE_FILE)
        cat("  check against 08: max | |CORE| - |CORE_08| | = ",
            format(max(abs(abs(fit$CORE) - abs(C08))), digits = 3), "\n", sep = "")
    }
    sel <- rev_select_core(fit$U_CELL, fit$CORE, fit$U_MODAL, meta$cell_type, groups)
    u <- rev_ep_factor_column(fit, S[vm], map, sel$a)
    gene <- rev_gene_scores(map, sel$G_oriented * u)
    gene[, variant := vn]

    e_ct <- rev_kw_epsilon2(sel$cell_loading_oriented, meta$cell_type)
    e_rep <- rev_kw_epsilon2(sel$cell_loading_oriented, meta$replicate)

    ord <- order(abs(as.numeric(fit$CORE)), decreasing = TRUE)[seq_len(min(TOP_CORE, length(fit$CORE)))]
    bs <- arrayInd(ord, dim(fit$CORE))[, 2]
    tb <- table(bs)
    top100 <- rbindlist(lapply(as.integer(names(tb)), function(b) data.table(
        variant = vn, CELL_component = sprintf("CELL%02d", b),
        n_top100_cores = as.integer(tb[as.character(b)]),
        celltype_epsilon2 = rev_kw_epsilon2(fit$U_CELL[, b], meta$cell_type)$epsilon2,
        replicate_epsilon2 = rev_kw_epsilon2(fit$U_CELL[, b], meta$replicate)$epsilon2)))

    go <- data.table()
    if (RUN_GSEA) {
        cat("  GSEA for ", vn, "...\n", sep = "")
        go <- rev_run_gsea(gene)
        if (nrow(go)) go[, variant := vn]
    }

    res[[vn]] <- list(
        gene = gene, go = go, top100 = top100,
        row = data.table(
            variant = vn, modalities = paste(vm, collapse = "+"),
            core_energy_fraction = fit$core_fro2 / fit$tensor_fro2,
            EP = sel$a, CELL = sel$b, K = sel$k, MOD = sel$r,
            core_rank_of_selected = sel$core_rank,
            core_value_oriented = sel$G_oriented,
            selection_score = sel$component_scores$score[1],
            abs_rank_biserial_A_vs_B = abs(sel$component_scores$rank_biserial_A_vs_B[1]),
            celltype_epsilon2 = e_ct$epsilon2, replicate_epsilon2 = e_rep$epsilon2,
            top100_n_components = nrow(top100),
            top100_n_celltype_gt_replicate = sum(top100$celltype_epsilon2 > top100$replicate_epsilon2),
            top100_median_celltype_eps2 = median(top100$celltype_epsilon2),
            top100_median_replicate_eps2 = median(top100$replicate_epsilon2)))
}

var_dt <- rbindlist(lapply(res, `[[`, "row"))
cat("\nVariants:\n")
print(var_dt[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)])

# ------------------------------------------------------------
# 2. Agreement with the four-modality result
# ------------------------------------------------------------

rank_agree <- rbindlist(lapply(setdiff(names(res), "ALL4"), function(vn)
    rev_compare_rankings(res$ALL4$gene, res[[vn]]$gene, "ALL4", vn, TOP_N)))
cat("\nGene-ranking agreement with ALL4:\n")
print(rank_agree[, .(variant = run_y, spearman = round(spearman, 4), side, top_n,
                     overlap, expected = round(expected, 2),
                     hypergeom_p = signif(hypergeom_p, 3))])

go_agree <- data.table(); go_top <- data.table()
if (RUN_GSEA) {
    go_agree <- rbindlist(lapply(setdiff(names(res), "ALL4"), function(vn)
        rev_compare_go(res$ALL4$go, res[[vn]]$go, "ALL4", vn, FDR)))
    go_top <- rbindlist(lapply(names(res), function(vn) rev_go_top_terms(res[[vn]]$go, vn, 10L)))
    cat("\nGO agreement with ALL4:\n")
    print(go_agree)
}

fwrite(var_dt, OUT_VAR, sep = "\t")
fwrite(rank_agree, OUT_RANK, sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "top100")), OUT_TOP100, sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "gene")), OUT_GENES, sep = "\t")
if (RUN_GSEA) {
    fwrite(go_agree, OUT_GO, sep = "\t")
    fwrite(go_top, OUT_GO_TOP, sep = "\t")
}
cat("\nSaved outputs with prefix GSE303006_M23_rev24_.\n")
rev_end_log(log_con)
