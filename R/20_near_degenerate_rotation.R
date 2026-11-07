# ============================================================
# 20_near_degenerate_rotation.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 1, Comment 2).
#
# The Stage-2 E-P and cell factors are invariant to any permutation
# of the Stage-1 component index k (this is exact and is shown
# algebraically in the Supplementary Note). Rotations inside a
# NEARLY degenerate singular subspace are not exactly neutral,
# because the singular values in the block differ slightly. This
# script measures how much such rotations change the results.
#
# Procedure
#   1. For each modality, find blocks of consecutive Stage-1
#      components whose adjacent singular values differ by less
#      than GAP_TOL (relative gap (d_k - d_{k+1}) / d_k).
#   2. Write the Stage-1 factors as  B = U diag(d~),  V.
#      For every block J draw a rotation Q (Haar-random orthogonal;
#      for 2-component blocks also the 45-degree rotation) and set
#         U_J <- U_J Q,   V_J <- V_J Q,   d~ unchanged.
#      This is the ambiguity of the singular vectors inside the
#      block; the reconstruction changes only through the small
#      differences of d~ inside the block.
#      All blocks of all modalities are rotated simultaneously.
#   3. Re-run the implicit Stage-2 HOSVD (revision_lib.R; same
#      algebra as 08) on the rotated Stage-1 factors.
#   4. Compare with the unrotated decomposition:
#        - largest principal angle between the retained cell
#          subspaces and between the retained E-P subspaces
#        - the core matched to the reference core (default: the
#          manuscript core #12; matched by maximal |correlation|
#          of the E-P and cell factor columns; k and r by maximal
#          |core|)
#        - Spearman correlation of the degree-corrected gene
#          scores and top-N gene overlap with the reference.
#
# Settings can be changed below or with environment variables
#   TDCHARM_GAP_TOL, TDCHARM_N_ROT
#
# Inputs (working directory): outputs of 06-09, 11
# Outputs:
#   GSE303006_M23_rev20_near_degenerate_blocks.tsv
#   GSE303006_M23_rev20_rotation_results.tsv
#   GSE303006_M23_rev20_rotation_summary.tsv
#   GSE303006_M23_rev20_rotation_log.txt
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

GAP_TOL <- as.numeric(Sys.getenv("TDCHARM_GAP_TOL", unset = "0.05"))
N_ROT <- as.integer(Sys.getenv("TDCHARM_N_ROT", unset = "20"))
REF_CORE_RANK <- 12L
TOP_N <- c(100L, 500L)
SEED <- 20261103L

MANIFEST_FILE <- "GSE303006_M23_stage2_implicit_manifest.rds"
MAPPING_FILE <- "GSE303006_M23_stage2_EP_mapping_degrees.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
EP_FACTOR_FILE <- "GSE303006_M23_stage2_Tucker_EP_factor.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"

OUT_BLOCKS <- "GSE303006_M23_rev20_near_degenerate_blocks.tsv"
OUT_RESULTS <- "GSE303006_M23_rev20_rotation_results.tsv"
OUT_SUMMARY <- "GSE303006_M23_rev20_rotation_summary.tsv"
OUT_LOG <- "GSE303006_M23_rev20_rotation_log.txt"

need <- c(MANIFEST_FILE, MAPPING_FILE, CORE_FILE, SMALL_FILE,
          EP_FACTOR_FILE, TOP_CORE_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "20_near_degenerate_rotation.R")
cat("GAP_TOL =", GAP_TOL, " N_ROT =", N_ROT, " SEED =", SEED, "\n\n")

# ------------------------------------------------------------
# 1. Load Stage-1 / Stage-2 objects
# ------------------------------------------------------------

mods <- c("ATAC", "H3K27ac", "RNA", "3D")
M <- readRDS(MANIFEST_FILE)
map <- as.data.table(readRDS(MAPPING_FILE))
S0 <- lapply(M$stage1_files, readRDS)[mods]

CORE08 <- readRDS(CORE_FILE)
SF08 <- readRDS(SMALL_FILE)
UEP08 <- readRDS(EP_FACTOR_FILE)

# Reference fit with the shared implementation; check it against 08
ref <- rev_stage2_hosvd(S0, map)
chk_core <- max(abs(abs(ref$CORE) - abs(CORE08)))
chk_cell <- max(abs(abs(crossprod(ref$U_CELL, SF08$U_CELL)) - diag(ncol(ref$U_CELL))))
cat("Check against 08: max | |CORE| - |CORE_08| | = ", format(chk_core, digits = 3),
    "; cell-factor alignment error = ", format(chk_cell, digits = 3), "\n", sep = "")
if (chk_core > 1e-8 || chk_cell > 1e-8) warning("Shared Stage-2 implementation differs from 08.")

top <- fread(TOP_CORE_FILE)
rc <- top[core_rank == REF_CORE_RANK]
ref_a <- rc$a_EP; ref_b <- rc$b_CELL; ref_k <- rc$k; ref_r <- rc$r_MODAL
cat(sprintf("Reference core rank %d: (EP%02d, CELL%02d, K%02d, MOD%02d)\n\n",
            REF_CORE_RANK, ref_a, ref_b, ref_k, ref_r))

ref_u <- UEP08[, ref_a]
ref_gene <- rev_gene_scores(map, CORE08[ref_a, ref_b, ref_k, ref_r] * ref_u)

# ------------------------------------------------------------
# 2. Near-degenerate blocks
# ------------------------------------------------------------

blocks <- list()
block_rows <- list()
for (m in mods) {
    d <- as.numeric(S0[[m]]$singular_values)
    K <- length(d)
    gap <- (d[-K] - d[-1]) / d[-K]
    link <- gap < GAP_TOL
    k <- 1L
    while (k <= K) {
        j <- k
        while (j < K && link[j]) j <- j + 1L
        if (j > k) {
            blocks[[length(blocks) + 1L]] <- list(modality = m, ks = k:j)
            block_rows[[length(block_rows) + 1L]] <- data.table(
                modality = m, k_from = k, k_to = j, size = j - k + 1L,
                d_first = d[k], d_last = d[j],
                max_rel_gap_in_block = max(gap[k:(j - 1L)]))
        }
        k <- j + 1L
    }
}
block_dt <- rbindlist(block_rows)
gap_dt <- rbindlist(lapply(mods, function(m) {
    d <- as.numeric(S0[[m]]$singular_values)
    data.table(modality = m, k = seq_len(length(d) - 1L),
               rel_gap_to_next = (d[-length(d)] - d[-1]) / d[-length(d)])
}))

cat("Relative gaps between consecutive Stage-1 singular values:\n")
print(dcast(gap_dt, k ~ modality, value.var = "rel_gap_to_next")[, lapply(.SD, signif, 3)])
cat("\nNear-degenerate blocks (relative gap < ", GAP_TOL, "):\n", sep = "")
if (nrow(block_dt)) print(block_dt) else cat("none\n")
fwrite(if (nrow(block_dt)) block_dt else data.table(note = "no near-degenerate blocks"),
       OUT_BLOCKS, sep = "\t")

# ------------------------------------------------------------
# 3. Rotation experiments
# ------------------------------------------------------------

rotate_stage1 <- function(S, rot) {
    for (bl in rot) {
        obj <- S[[bl$modality]]
        d_scaled <- obj$singular_values / sqrt(sum(obj$singular_values^2))
        B <- obj$feature_scores
        U <- sweep(B, 2, d_scaled, "/")
        V <- obj$cell_loadings
        J <- bl$ks
        U[, J] <- U[, J, drop = FALSE] %*% bl$Q
        V[, J] <- V[, J, drop = FALSE] %*% bl$Q
        obj$feature_scores <- sweep(U, 2, d_scaled, "*")
        obj$cell_loadings <- V
        S[[bl$modality]] <- obj
    }
    S
}

haar <- function(n) {
    z <- matrix(rnorm(n * n), n, n)
    q <- qr(z); Q <- qr.Q(q)
    Q %*% diag(sign(diag(qr.R(q))), n)
}

max_angle_deg <- function(A, B) {
    sv <- svd(crossprod(A, B))$d
    acos(min(1, max(0, min(sv)))) * 180 / pi
}

designs <- list()
if (length(blocks)) {
    set.seed(SEED)
    for (i in seq_len(N_ROT)) {
        designs[[length(designs) + 1L]] <- list(
            id = sprintf("haar_%02d", i),
            rot = lapply(blocks, function(bl) c(bl, list(Q = haar(length(bl$ks))))))
    }
    if (any(vapply(blocks, function(bl) length(bl$ks) == 2L, logical(1)))) {
        th <- pi / 4
        designs[[length(designs) + 1L]] <- list(
            id = "deg45_all_pairs",
            rot = lapply(Filter(function(bl) length(bl$ks) == 2L, blocks), function(bl)
                c(bl, list(Q = matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)))))
    }
}

res_rows <- list()

for (ds in designs) {

    S1 <- rotate_stage1(S0, ds$rot)
    fit <- rev_stage2_hosvd(S1, map, verbose = FALSE)

    # cell subspace angle
    ang_cell <- max_angle_deg(ref$U_CELL, fit$U_CELL)

    # E-P subspace angle: U = A W Sigma^-1 for each fit.
    # U_ref^T U_fit = Sigma_ref^-1 W_ref^T A_ref^T A_fit W_fit Sigma_fit^-1
    P <- nrow(map)
    bl <- split(seq_len(P), ceiling(seq_len(P) / 10000L))
    HA_x <- matrix(0, nrow(ref$H_A), nrow(fit$H_A))
    for (ii in bl) {
        HA_x <- HA_x + crossprod(rev_build_A_block(S0, map, ii),
                                 rev_build_A_block(S1, map, ii))
    }
    cross <- diag(1 / ref$sigmaP) %*% t(ref$WP) %*% HA_x %*% fit$WP %*% diag(1 / fit$sigmaP)
    ang_ep <- acos(min(1, max(0, min(svd(cross)$d)))) * 180 / pi

    # matched core
    a <- which.max(abs(cross[ref_a, ]))
    b <- which.max(abs(cor(ref$U_CELL[, ref_b], fit$U_CELL)))
    sub <- fit$CORE[a, b, , , drop = FALSE]
    idx <- arrayInd(which.max(abs(as.numeric(sub))), dim(sub))
    k <- idx[3]; r <- idx[4]
    u <- rev_ep_factor_column(fit, S1, map, a)
    gene <- rev_gene_scores(map, fit$CORE[a, b, k, r] * u)
    cmp <- rev_compare_rankings(ref_gene, gene, "reference", ds$id, TOP_N)
    # gene-score sign is arbitrary only through matched-column signs;
    # report |rho| and orient overlaps accordingly
    sgn <- sign(cmp$spearman[1])
    if (sgn < 0) {
        gene[, gene_score := -gene_score]
        cmp <- rev_compare_rankings(ref_gene, gene, "reference", ds$id, TOP_N)
    }

    ov <- dcast(cmp, run_y ~ paste0("top", top_n, "_overlap_",
                                    ifelse(side == "positive", "pos", "neg")),
                value.var = "overlap")[, run_y := NULL]
    res_rows[[length(res_rows) + 1L]] <- cbind(data.table(
        design = ds$id,
        max_angle_cell_subspace_deg = ang_cell,
        max_angle_EP_subspace_deg = ang_ep,
        matched_EP = a, matched_CELL = b, matched_K = k, matched_MOD = r,
        abs_corr_EP = abs(cross[ref_a, a]),
        abs_corr_CELL = abs(cor(ref$U_CELL[, ref_b], fit$U_CELL[, b])),
        gene_spearman = abs(cmp$spearman[1])), ov)
    cat(sprintf("%-16s cell angle %.3f deg, EP angle %.3f deg, matched core (%d,%d,%d,%d), gene rho %.4f\n",
                ds$id, ang_cell, ang_ep, a, b, k, r, abs(cmp$spearman[1])))
}

if (length(res_rows)) {
    res <- rbindlist(res_rows)
    ovc <- grep("_overlap_", names(res), value = TRUE)
    summ <- cbind(
        res[, .(n_designs = .N,
                max_angle_cell_deg_max = max(max_angle_cell_subspace_deg),
                max_angle_EP_deg_max = max(max_angle_EP_subspace_deg),
                abs_corr_EP_min = min(abs_corr_EP),
                abs_corr_CELL_min = min(abs_corr_CELL),
                gene_spearman_min = min(gene_spearman),
                gene_spearman_median = median(gene_spearman),
                same_K_as_reference = sum(matched_K == ref_k))],
        res[, lapply(.SD, min), .SDcols = ovc][, setnames(.SD, paste0(names(.SD), "_min"))])
    cat("\nSummary over rotation designs:\n")
    print(t(summ))
} else {
    res <- data.table(note = "no near-degenerate blocks; no rotation performed")
    summ <- res
    cat("\nNo near-degenerate blocks at GAP_TOL = ", GAP_TOL, "; nothing to rotate.\n", sep = "")
}

fwrite(res, OUT_RESULTS, sep = "\t")
fwrite(summ, OUT_SUMMARY, sep = "\t")
cat("\nSaved:\n", OUT_BLOCKS, "\n", OUT_RESULTS, "\n", OUT_SUMMARY, "\n", OUT_LOG, "\n", sep = "")
rev_end_log(log_con)
