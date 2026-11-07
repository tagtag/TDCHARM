# ============================================================
# 22_celltype_replicate_joint.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 1, Comment 4).
#
# The manuscript compares cell-type and replicate effects on the
# cell components marginally (two separate Kruskal-Wallis tests).
# This script provides
#
#   1. the cell-type x replicate composition of the 4,258 cells
#      (counts, row proportions, chi-square test of independence);
#   2. for every cell component (all 20; flagged if it appears in the
#      top-100 core elements) a joint linear model
#
#        loading ~ cell_type + replicate + cell_type:replicate
#                  + log(ATAC fragments) + log(H3K27ac fragments)
#                  + log(RNA library size) + log(RNA detected genes)
#                  + 3D observed fraction
#
#      with type-II sums of squares and partial eta^2 for every term
#      (partial eta^2 = SS_term / (SS_term + SS_residual of the full
#      model)). The same model is also fitted to rank-transformed
#      loadings as a robustness check;
#   3. within-cell-type replicate effects: for every component and
#      every cell type with at least MIN_CELLS_PER_REP cells in at
#      least two replicates, a Kruskal-Wallis test of the loading
#      across replicates (epsilon^2, BH-adjusted P over all tests).
#
# The 3D observed fraction of a cell is the fraction of the final
# 20-kb bin pairs whose reconstructed distance is observed (finite
# and positive) in that cell. It is computed once from the 3D HDF5
# file and cached.
#
# Inputs (working directory): outputs of 03-05, 08, 09, 10
# Outputs:
#   GSE303006_M23_rev22_celltype_by_replicate_counts.tsv
#   GSE303006_M23_rev22_celltype_by_replicate_proportions.tsv
#   GSE303006_M23_rev22_cell_covariates.tsv.gz
#   GSE303006_M23_rev22_joint_model_partial_eta2.tsv
#   GSE303006_M23_rev22_joint_model_summary.tsv
#   GSE303006_M23_rev22_within_celltype_replicate_KW.tsv
#   GSE303006_M23_rev22_within_celltype_replicate_summary.tsv
#   GSE303006_M23_rev22_log.txt
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

MIN_CELLS_PER_REP <- 20L
TOP_CORE <- 100L
BLOCK_3D <- 2000L

META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"
QC_FILE <- "GSE303006_M23_brain_final_cell_QC.tsv"
EP_FILE <- "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"
H5_3D <- "GSE303006_M23_brain_3D_binpair_by_cell.h5"
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"

OUT_COUNTS <- "GSE303006_M23_rev22_celltype_by_replicate_counts.tsv"
OUT_PROP <- "GSE303006_M23_rev22_celltype_by_replicate_proportions.tsv"
OUT_COV <- "GSE303006_M23_rev22_cell_covariates.tsv.gz"
OUT_ETA <- "GSE303006_M23_rev22_joint_model_partial_eta2.tsv"
OUT_ETA_SUM <- "GSE303006_M23_rev22_joint_model_summary.tsv"
OUT_WITHIN <- "GSE303006_M23_rev22_within_celltype_replicate_KW.tsv"
OUT_WITHIN_SUM <- "GSE303006_M23_rev22_within_celltype_replicate_summary.tsv"
OUT_LOG <- "GSE303006_M23_rev22_log.txt"

need <- c(META_FILE, QC_FILE, EP_FILE, SMALL_FILE, TOP_CORE_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "22_celltype_replicate_joint.R")

SF <- readRDS(SMALL_FILE)
U <- SF$U_CELL
cells <- rownames(U)
meta <- fread(META_FILE)[match(cells, cell)]
stopifnot(identical(meta$cell, cells))

top <- fread(TOP_CORE_FILE)
top_b <- sort(unique(top[core_rank <= TOP_CORE, b_CELL]))

# ------------------------------------------------------------
# 1. Composition
# ------------------------------------------------------------

tab <- table(cell_type = meta$cell_type, replicate = meta$replicate, useNA = "no")
counts <- as.data.table(as.data.frame.matrix(tab), keep.rownames = "cell_type")
counts[, total := rowSums(.SD), .SDcols = setdiff(names(counts), "cell_type")]
setorder(counts, -total)
prop <- copy(counts)
reps <- setdiff(names(counts), c("cell_type", "total"))
prop[, (reps) := lapply(.SD, function(x) x / total), .SDcols = reps]
chi <- suppressWarnings(chisq.test(tab))
cramer <- sqrt(as.numeric(chi$statistic) / (sum(tab) * (min(dim(tab)) - 1)))

cat("Cell-type x replicate composition (counts):\n")
print(counts)
cat(sprintf("\nChi-square test of independence: X2 = %.2f, df = %d, P = %.3g, Cramer's V = %.3f\n",
            chi$statistic, chi$parameter, chi$p.value, cramer))

# ------------------------------------------------------------
# 2. Covariates (depth / quality)
# ------------------------------------------------------------

qc <- fread(QC_FILE)[match(cells, cell)]
cov <- data.table(cell = cells,
                  log_ATAC = log(qc$ATAC_count),
                  log_H3K27ac = log(qc$H3K27ac_count),
                  log_RNA_lib = log(qc$RNA_library_size),
                  log_RNA_genes = log(qc$RNA_detected_genes))

if (file.exists(OUT_COV)) {
    old <- fread(OUT_COV)
    if (identical(old$cell, cells) && "frac_3D_observed" %in% names(old)) {
        cov[, frac_3D_observed := old$frac_3D_observed]
        cat("\n3D observed fraction read from cache.\n")
    }
}
if (!"frac_3D_observed" %in% names(cov)) {
    if (!file.exists(H5_3D) || !requireNamespace("rhdf5", quietly = TRUE)) {
        warning("3D HDF5 or rhdf5 unavailable; 3D covariate omitted.")
    } else {
        cat("\nComputing 3D observed fraction per cell from the HDF5 file...\n")
        bp <- sort(unique(as.integer(readRDS(EP_FILE)$binpair_id)))
        h5_bp <- as.integer(rhdf5::h5read(H5_3D, "binpair_id"))
        h5_cells <- as.character(rhdf5::h5read(H5_3D, "cell"))
        rows <- match(bp, h5_bp); cols <- match(cells, h5_cells)
        stopifnot(!anyNA(rows), !anyNA(cols))
        obs <- numeric(length(cells))
        for (ii in split(seq_along(rows), ceiling(seq_along(rows) / BLOCK_3D))) {
            d <- rhdf5::h5read(H5_3D, "distance", index = list(rows[ii], cols))
            obs <- obs + colSums(is.finite(d) & d > 0)
        }
        rhdf5::h5closeAll()
        cov[, frac_3D_observed := obs / length(rows)]
    }
}
fwrite(cov, OUT_COV, sep = "\t")
cat("\nCovariate summary:\n")
print(summary(cov[, -1]))

covars <- intersect(c("log_ATAC", "log_H3K27ac", "log_RNA_lib", "log_RNA_genes",
                      "frac_3D_observed"), names(cov))

# ------------------------------------------------------------
# 3. Joint model with type-II SS
# ------------------------------------------------------------

fml <- function(...) {
    rhs <- unlist(list(...))
    rhs <- rhs[nzchar(rhs)]
    as.formula(paste("y ~", paste(rhs, collapse = " + ")))
}

type2_partial_eta2 <- function(y, d) {
    d$y <- y
    f_full <- fml("cell_type * replicate", covars)
    f_main <- fml("cell_type + replicate", covars)
    fit_full <- lm(f_full, data = d)
    rss_full <- sum(residuals(fit_full)^2)
    df_res <- fit_full$df.residual
    # type II: each main effect adjusted for the other main effect and
    # the covariates; the interaction adjusted for both main effects;
    # each covariate adjusted for everything else.
    terms <- list(
        cell_type = list(fml("replicate", covars), f_main),
        replicate = list(fml("cell_type", covars), f_main),
        `cell_type:replicate` = list(f_main, f_full))
    for (v in covars) {
        terms[[v]] <- list(fml("cell_type * replicate", setdiff(covars, v)), f_full)
    }
    rbindlist(lapply(names(terms), function(tn) {
        f0 <- terms[[tn]][[1]]; f1 <- terms[[tn]][[2]]
        m0 <- lm(f0, data = d); m1 <- lm(f1, data = d)
        ss <- sum(residuals(m0)^2) - sum(residuals(m1)^2)
        df <- m0$df.residual - m1$df.residual
        Fv <- (ss / df) / (rss_full / df_res)
        data.table(term = tn, df = df, SS = ss,
                   partial_eta2 = ss / (ss + rss_full),
                   F = Fv, p = pf(Fv, df, df_res, lower.tail = FALSE))
    }))
}

d0 <- data.table(cell_type = factor(meta$cell_type), replicate = factor(meta$replicate))
d0 <- cbind(d0, cov[, ..covars])
ok <- complete.cases(d0)
cat(sprintf("\nCells used in the joint model: %d / %d\n", sum(ok), nrow(d0)))
d <- as.data.frame(d0[ok])

eta <- rbindlist(lapply(seq_len(ncol(U)), function(b) {
    y <- U[ok, b]
    rbind(
        cbind(data.table(CELL_component = colnames(U)[b], response = "loading"),
              type2_partial_eta2(y, d)),
        cbind(data.table(CELL_component = colnames(U)[b], response = "rank(loading)"),
              type2_partial_eta2(rank(y), d)))
}))
eta[, in_top100 := as.integer(sub("CELL", "", CELL_component)) %in% top_b]

summ <- dcast(eta, CELL_component + in_top100 + response ~ term,
              value.var = "partial_eta2")
setcolorder(summ, c("CELL_component", "in_top100", "response", "cell_type",
                    "replicate", "cell_type:replicate"))
summ[, celltype_gt_replicate := cell_type > replicate]
summ[, celltype_gt_interaction := cell_type > `cell_type:replicate`]

cat("\nPartial eta^2 (type-II SS), raw loadings, top-100 components:\n")
print(summ[in_top100 == TRUE & response == "loading",
           lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat("\nCounts over top-100 components (raw / rank):\n")
print(summ[in_top100 == TRUE, .(n = .N,
                                celltype_gt_replicate = sum(celltype_gt_replicate),
                                celltype_gt_interaction = sum(celltype_gt_interaction),
                                median_eta2_celltype = median(cell_type),
                                median_eta2_replicate = median(replicate),
                                median_eta2_interaction = median(`cell_type:replicate`)),
           by = response])

# ------------------------------------------------------------
# 4. Within-cell-type replicate effects
# ------------------------------------------------------------

ct_ok <- counts[, c("cell_type", reps), with = FALSE]
ct_ok[, n_reps_ok := rowSums(.SD >= MIN_CELLS_PER_REP), .SDcols = reps]
eligible <- ct_ok[n_reps_ok >= 2L, cell_type]
cat(sprintf("\nCell types with >= %d cells in >= 2 replicates: %d of %d\n",
            MIN_CELLS_PER_REP, length(eligible), nrow(ct_ok)))

within <- rbindlist(lapply(seq_len(ncol(U)), function(b) {
    rbindlist(lapply(eligible, function(ct) {
        idx <- which(meta$cell_type == ct)
        rr <- meta$replicate[idx]
        keep_r <- names(which(table(rr) >= MIN_CELLS_PER_REP))
        idx <- idx[rr %in% keep_r]
        e <- rev_kw_epsilon2(U[idx, b], meta$replicate[idx])
        data.table(CELL_component = colnames(U)[b], cell_type = ct,
                   replicates_tested = paste(sort(keep_r), collapse = ","),
                   n_cells = length(idx), H = e$statistic, df = e$df,
                   p = e$p, epsilon2 = e$epsilon2)
    }))
}))
within[, p_BH := p.adjust(p, "BH")]
within[, in_top100 := as.integer(sub("CELL", "", CELL_component)) %in% top_b]

wsum <- within[, .(n_tests = .N, n_BH_005 = sum(p_BH <= 0.05),
                   median_epsilon2 = median(epsilon2),
                   max_epsilon2 = max(epsilon2),
                   cell_type_at_max = cell_type[which.max(epsilon2)]),
               by = .(CELL_component, in_top100)]
cat("\nWithin-cell-type replicate Kruskal-Wallis, per component:\n")
print(wsum[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------

fwrite(counts, OUT_COUNTS, sep = "\t")
fwrite(prop, OUT_PROP, sep = "\t")
fwrite(eta, OUT_ETA, sep = "\t")
fwrite(summ, OUT_ETA_SUM, sep = "\t")
fwrite(within, OUT_WITHIN, sep = "\t")
fwrite(wsum, OUT_WITHIN_SUM, sep = "\t")
cat("\nSaved outputs with prefix GSE303006_M23_rev22_.\n")
rev_end_log(log_con)
