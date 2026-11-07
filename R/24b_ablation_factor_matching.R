# ============================================================
# 24b_ablation_factor_matching.R
# GSE303006 / GENCODE M23
#
# Follow-up to 24_modality_ablation.R (Reviewer 1, Comment 1), in the
# same way as 21b/21c follow up 21.
#
# In 24 the pre-specified selection rule picked, in some variants, a
# core element far down the core ranking (e.g. rank 116 without ATAC).
# That alone does not tell whether the E-P program of core #12 is
# absent from the variant or merely not picked by the rule. Here the
# factor VECTORS of every variant are compared directly with those of
# the four-modality decomposition, without any selection rule:
#
#   E-P factors  : |Pearson r| over the 730,969 E-P pairs. Because
#                  U_EP = A W Sigma^{-1} with A the degree-corrected
#                  Stage-1 feature scores (as in revision_lib.R), the
#                  cross-products are accumulated block-wise in one
#                  pass; no factor matrix is stored.
#   cell factors : |Pearson r| over the 4,258 cells
#
# For each variant, the factors best matching the E-P factor (a) and
# the cell factor (b) of core #12 are taken, and the core element with
# these two factors and the largest |core| over (k, modality) is
# analysed as in 21c: core value and rank, degree-corrected gene
# ranking and its agreement with core #12, and GO-BP GSEA (as in 12;
# serial; TDCHARM_RUN_GSEA=0 skips it).
#
# Orientation. The gene score G u_EP is invariant to the sign of the
# E-P factor; it depends on the signs of the cell and modality
# columns. As in 21, the cell column is oriented so that group A is
# higher. The modality column is oriented by its inner product with
# the core-#12 modality column of the four-modality fit (itself
# oriented to a positive sum, as in 21) over the modalities the variant
# contains; with an unchanged modality set this equals the rule of 21,
# whereas a positive-sum rule is not comparable between different
# modality sets.
#
# Inputs (working directory): outputs of 06-10, and
#   GSE303006_M23_rev24_variants.tsv (from 24; optional, for the
#   comparison with the rule-selected core)
# Outputs:
#   GSE303006_M23_rev24b_EP_factor_abs_correlation.tsv
#   GSE303006_M23_rev24b_CELL_factor_abs_correlation.tsv
#   GSE303006_M23_rev24b_best_matches.tsv
#   GSE303006_M23_rev24b_core12_matched.tsv
#   GSE303006_M23_rev24b_gene_ranking_agreement.tsv
#   GSE303006_M23_rev24b_gene_scores.tsv.gz
#   GSE303006_M23_rev24b_GO_agreement.tsv
#   GSE303006_M23_rev24b_GO_top10_per_variant.tsv
#   GSE303006_M23_rev24b_log.txt
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

CORE_RANK <- 12L
TOP_N <- c(100L, 500L)
FDR <- 0.05
PAIR_BLOCK <- 10000L
RUN_GSEA <- !identical(Sys.getenv("TDCHARM_RUN_GSEA"), "0") &&
    requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Mm.eg.db", quietly = TRUE)

MANIFEST_FILE <- "GSE303006_M23_stage2_implicit_manifest.rds"
MAPPING_FILE <- "GSE303006_M23_stage2_EP_mapping_degrees.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"
VAR24_FILE <- "GSE303006_M23_rev24_variants.tsv"

OUT_EP <- "GSE303006_M23_rev24b_EP_factor_abs_correlation.tsv"
OUT_CELL <- "GSE303006_M23_rev24b_CELL_factor_abs_correlation.tsv"
OUT_BEST <- "GSE303006_M23_rev24b_best_matches.tsv"
OUT_CORE <- "GSE303006_M23_rev24b_core12_matched.tsv"
OUT_RANK <- "GSE303006_M23_rev24b_gene_ranking_agreement.tsv"
OUT_GENES <- "GSE303006_M23_rev24b_gene_scores.tsv.gz"
OUT_GO <- "GSE303006_M23_rev24b_GO_agreement.tsv"
OUT_GO_TOP <- "GSE303006_M23_rev24b_GO_top10_per_variant.tsv"
OUT_LOG <- "GSE303006_M23_rev24b_log.txt"

need <- c(MANIFEST_FILE, MAPPING_FILE, CORE_FILE, TOP_CORE_FILE, META_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "24b_ablation_factor_matching.R")
cat("GSEA: ", if (RUN_GSEA) "run" else "skipped", "\n\n", sep = "")

mods <- c("ATAC", "H3K27ac", "RNA", "3D")
M <- readRDS(MANIFEST_FILE)
map <- as.data.table(readRDS(MAPPING_FILE))
S <- lapply(M$stage1_files, readRDS)[mods]
cells <- rownames(S[[1]]$cell_loadings)
meta <- fread(META_FILE)[match(cells, cell)]
groups <- rev_contrast_groups(meta$cell_type)
inA <- meta$cell_type %in% groups$A
inB <- meta$cell_type %in% groups$B

variants <- c(list(ALL4 = mods),
              setNames(lapply(mods, function(m) setdiff(mods, m)), paste0("minus_", mods)),
              setNames(as.list(mods), paste0("only_", mods)))

top <- fread(TOP_CORE_FILE)
cr <- top[core_rank == CORE_RANK]
a0 <- cr$a_EP; b0 <- cr$b_CELL; k0 <- cr$k; r0 <- cr$r_MODAL
cat(sprintf("Reference: four-modality core #%d = (EP%02d, CELL%02d, K%02d, MOD%02d)\n\n",
            CORE_RANK, a0, b0, k0, r0))

# ------------------------------------------------------------
# 1. Stage 2 for every variant (as in 24)
# ------------------------------------------------------------

fits <- lapply(variants, function(vm) rev_stage2_hosvd(S[vm], map))
C08 <- readRDS(CORE_FILE)
cat("Check ALL4 against 08: max | |CORE| - |CORE_08| | = ",
    format(max(abs(abs(fits$ALL4$CORE) - abs(C08))), digits = 3), "\n\n", sep = "")
rm(C08)

# Columns of the full feature matrix A (all four modalities) used by each variant
K <- fits$ALL4$K
off_full <- setNames(lapply(seq_along(mods), function(i) ((i - 1L) * K + 1L):(i * K)), mods)
cols <- lapply(variants, function(vm) unlist(off_full[vm], use.names = FALSE))
proj <- lapply(fits, function(f) sweep(f$WP, 2, f$sigmaP, "/"))   # L_v x R_EP

# ------------------------------------------------------------
# 2. E-P factor correlations, accumulated over pair blocks
# ------------------------------------------------------------

P <- nrow(map)
vn_all <- names(variants)
sum_u <- lapply(fits, function(f) numeric(ncol(f$WP)))
ssq_u <- sum_u
xprod <- lapply(fits, function(f) matrix(0, ncol(fits$ALL4$WP), ncol(f$WP)))
for (ii in split(seq_len(P), ceiling(seq_len(P) / PAIR_BLOCK))) {
    A <- rev_build_A_block(S, map, ii)              # columns in the order of mods
    U0 <- A %*% proj$ALL4
    for (vn in vn_all) {
        Uv <- if (vn == "ALL4") U0 else A[, cols[[vn]], drop = FALSE] %*% proj[[vn]]
        sum_u[[vn]] <- sum_u[[vn]] + colSums(Uv)
        ssq_u[[vn]] <- ssq_u[[vn]] + colSums(Uv^2)
        xprod[[vn]] <- xprod[[vn]] + crossprod(U0, Uv)
    }
}
pearson <- function(vn) {
    m0 <- sum_u$ALL4 / P; mv <- sum_u[[vn]] / P
    cv <- xprod[[vn]] - P * tcrossprod(m0, mv)
    s0 <- sqrt(ssq_u$ALL4 - P * m0^2); sv <- sqrt(ssq_u[[vn]] - P * mv^2)
    R <- abs(cv / tcrossprod(s0, sv))
    dimnames(R) <- list(sprintf("EP%02d", seq_len(nrow(R))), sprintf("EP%02d", seq_len(ncol(R))))
    R
}
R_ep <- lapply(setNames(vn_all, vn_all), pearson)
cat("Self-check: |r| of ALL4 E-P factors with themselves (diagonal min) = ",
    format(min(diag(R_ep$ALL4)), digits = 6), "\n", sep = "")

# ------------------------------------------------------------
# 3. Cell factor correlations
# ------------------------------------------------------------

R_cell <- lapply(fits, function(f) {
    R <- abs(cor(fits$ALL4$U_CELL, f$U_CELL))
    dimnames(R) <- list(colnames(fits$ALL4$U_CELL), colnames(f$U_CELL)); R
})

best_table <- function(R, mode, vn) {
    rbindlist(lapply(seq_len(nrow(R)), function(i) {
        o <- order(R[i, ], decreasing = TRUE)
        data.table(mode = mode, variant = vn, ALL4_factor = rownames(R)[i],
                   best_factor = colnames(R)[o[1]], best_abs_r = R[i, o[1]],
                   second_factor = colnames(R)[o[2]], second_abs_r = R[i, o[2]])
    }))
}
best <- rbindlist(lapply(vn_all, function(vn)
    rbind(best_table(R_ep[[vn]], "EP", vn), best_table(R_cell[[vn]], "CELL", vn))))

ep0 <- sprintf("EP%02d", a0); cell0 <- sprintf("CELL%02d", b0)
cat("\nBest matches of the core-#12 factors in every variant:\n")
print(best[(mode == "EP" & ALL4_factor == ep0) | (mode == "CELL" & ALL4_factor == cell0),
           .(mode, variant, ALL4_factor, best_factor, best_abs_r = round(best_abs_r, 3),
             second_factor, second_abs_r = round(second_abs_r, 3))])

# ------------------------------------------------------------
# 4. Core element built from the matched factors, gene ranking, GO
# ------------------------------------------------------------

var24 <- if (file.exists(VAR24_FILE)) fread(VAR24_FILE) else NULL
w0 <- fits$ALL4$U_MODAL[, r0]
w0 <- w0 * (if (sum(w0) >= 0) 1 else -1)

res <- list()
for (vn in vn_all) {
    f <- fits[[vn]]; vm <- variants[[vn]]
    a <- as.integer(sub("EP", "", best[mode == "EP" & variant == vn & ALL4_factor == ep0, best_factor]))
    b <- as.integer(sub("CELL", "", best[mode == "CELL" & variant == vn & ALL4_factor == cell0, best_factor]))
    slice <- f$CORE[a, b, , , drop = FALSE]
    w <- arrayInd(which.max(abs(as.numeric(slice))), dim(slice))
    k <- w[3]; r <- w[4]
    G <- f$CORE[a, b, k, r]
    s_b <- if (mean(f$U_CELL[inA, b]) >= mean(f$U_CELL[inB, b])) 1 else -1
    s_r <- if (sum(f$U_MODAL[vm, r] * w0[vm]) >= 0) 1 else -1
    G_or <- G * s_b * s_r
    u <- rev_ep_factor_column(f, S[vm], map, a)
    gene <- rev_gene_scores(map, G_or * u)
    gene[, variant := vn]
    slice_b <- abs(as.numeric(f$CORE[, b, , ]))
    e_ct <- rev_kw_epsilon2(f$U_CELL[, b] * s_b, meta$cell_type)
    rb <- rev_rank_biserial(f$U_CELL[inA, b] * s_b, f$U_CELL[inB, b] * s_b)

    sel24 <- if (!is.null(var24)) var24[variant == vn] else NULL
    go <- data.table()
    if (RUN_GSEA) {
        cat("  GSEA for ", vn, "...\n", sep = "")
        go <- rev_run_gsea(gene)
        if (nrow(go)) go[, variant := vn]
    }
    res[[vn]] <- list(gene = gene, go = go, row = data.table(
        variant = vn, modalities = paste(vm, collapse = "+"),
        EP = a, EP_abs_r = R_ep[[vn]][ep0, a],
        CELL = b, CELL_abs_r = R_cell[[vn]][cell0, b],
        K = k, MOD = r, core_value_oriented = G_or,
        MOD_abs_cos_with_core12 = abs(sum(f$U_MODAL[vm, r] * w0[vm])) / sqrt(sum(w0[vm]^2)),
        rank_among_all_core = sum(abs(as.numeric(f$CORE)) > abs(G)) + 1L,
        rank_among_CELL_elements = sum(slice_b > abs(G)) + 1L,
        ratio_to_largest_core = abs(G) / max(abs(f$CORE)),
        rank_biserial_A_vs_B = rb, celltype_epsilon2 = e_ct$epsilon2,
        rule24_EP = if (length(sel24$EP)) sel24$EP else NA_real_,
        rule24_CELL = if (length(sel24$CELL)) sel24$CELL else NA_real_,
        rule24_core_rank = if (length(sel24$core_rank_of_selected)) sel24$core_rank_of_selected else NA_real_))
}

core_dt <- rbindlist(lapply(res, `[[`, "row"))
cat("\nCore element from the factors matching core #12, per variant\n",
    "(rule24_*: the element selected by the rule in 24, for comparison):\n", sep = "")
print(core_dt[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)])

rank_agree <- rbindlist(lapply(setdiff(vn_all, "ALL4"), function(vn)
    rev_compare_rankings(res$ALL4$gene, res[[vn]]$gene, "ALL4", vn, TOP_N)))
cat("\nGene-ranking agreement with core #12 (ALL4):\n")
print(rank_agree[, .(variant = run_y, spearman = round(spearman, 4), side, top_n,
                     overlap, expected = round(expected, 2),
                     hypergeom_p = signif(hypergeom_p, 3))])

go_agree <- data.table(); go_top <- data.table()
if (RUN_GSEA) {
    go_agree <- rbindlist(lapply(setdiff(vn_all, "ALL4"), function(vn)
        rev_compare_go(res$ALL4$go, res[[vn]]$go, "ALL4", vn, FDR)))
    go_top <- rbindlist(lapply(vn_all, function(vn) rev_go_top_terms(res[[vn]]$go, vn, 10L)))
    cat("\nGO agreement with core #12 (ALL4):\n")
    print(go_agree)
}

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------

fwrite(rbindlist(lapply(vn_all, function(vn)
    cbind(variant = vn, as.data.table(R_ep[[vn]], keep.rownames = "ALL4_factor"))), fill = TRUE), OUT_EP, sep = "\t")
fwrite(rbindlist(lapply(vn_all, function(vn)
    cbind(variant = vn, as.data.table(R_cell[[vn]], keep.rownames = "ALL4_factor"))), fill = TRUE), OUT_CELL, sep = "\t")
fwrite(best, OUT_BEST, sep = "\t")
fwrite(core_dt, OUT_CORE, sep = "\t")
fwrite(rank_agree, OUT_RANK, sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "gene")), OUT_GENES, sep = "\t")
if (RUN_GSEA) {
    fwrite(go_agree, OUT_GO, sep = "\t")
    fwrite(go_top, OUT_GO_TOP, sep = "\t")
}
cat("\nSaved outputs with prefix GSE303006_M23_rev24b_.\n")
rev_end_log(log_con)
