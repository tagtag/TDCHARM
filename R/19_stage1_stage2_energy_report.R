# ============================================================
# 19_stage1_stage2_energy_report.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 1, Comment 6; Reviewer 3, Comment 1):
#   Stage-1 variance retention for each modality and
#   reconstruction errors at both stages.
#
# No decomposition is recomputed. Everything is read from the
# saved outputs of 05 (Stage 1) and 08 (Stage 2).
#
# ------------------------------------------------------------
# Stage 1 (per modality m; Z_m = standardized feature x cell matrix)
#
#   ||Z_m||_F^2 = total energy
#       ATAC/H3K27ac/RNA: (#nonconstant features) x (n_cells - 1)
#       3D:               sum over nonconstant bin pairs of (n_obs - 1)
#     (identical to the denominator of explained_fraction in 05)
#
#   retained fraction  f1_m = sum_{k=1}^{K} d_mk^2 / ||Z_m||_F^2
#   relative error     e1_m = ||Z_m - Zhat_m||_F / ||Z_m||_F
#                           = sqrt(1 - f1_m)
#   (exact for truncated SVD; also exact for the randomized SVD of
#    3D, because Zhat = Q B_K is an orthogonal projection of Z.)
#
# Stage 2 (conceptual tensor X; Tucker factors orthonormal)
#
#   ||X - Xhat||_F^2 = ||X||_F^2 - ||G||_F^2
#   retained fraction  f2   = ||G||_F^2 / ||X||_F^2
#   relative error     e2   = sqrt(1 - f2)
#
#   Per-modality Stage-2 retention:
#   because the K mode (identity, 20/20) and the modality mode
#   (orthogonal, 4/4) are kept at full rank, rotating the core
#   back along the modality mode gives the projection of each
#   modality slice X_m onto the retained E-P and cell subspaces:
#
#     Gpre[,,,m] = sum_r CORE[,,,r] * U_MODAL[m, r]
#     f2_m = ||Gpre[,,,m]||_F^2 / ||X_m||_F^2
#
#   with sum_m ||Gpre[,,,m]||^2 = ||G||^2 (checked below).
#
# Overall (original standardized assay -> Stage-2 representation):
#     f_total_m = f1_m * f2_m
#   (the degree-corrected E-P mapping preserves Frobenius norms;
#    ||X_m||_F^2 is reported to verify this.)
#
# Outputs:
#   GSE303006_M23_rev_stage1_energy_by_modality.tsv
#   GSE303006_M23_rev_stage1_cumulative_energy.tsv
#   GSE303006_M23_rev_stage2_energy.tsv
#   GSE303006_M23_rev_stage2_energy_by_modality.tsv
#   GSE303006_M23_rev_stage1_stage2_energy_log.txt
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
# 0. Files
# ------------------------------------------------------------

STAGE1_FILES <- c(
    ATAC = "GSE303006_M23_stage1_K20_ATAC_signaligned.rds",
    H3K27ac = "GSE303006_M23_stage1_K20_H3K27ac_signaligned.rds",
    RNA = "GSE303006_M23_stage1_K20_RNA_signaligned.rds",
    `3D` = "GSE303006_M23_stage1_K20_3D_signaligned.rds"
)

DIAG_FILE <- "GSE303006_M23_stage2_Tucker_diagnostics.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"

OUT_S1 <- "GSE303006_M23_rev_stage1_energy_by_modality.tsv"
OUT_S1_CUM <- "GSE303006_M23_rev_stage1_cumulative_energy.tsv"
OUT_S2 <- "GSE303006_M23_rev_stage2_energy.tsv"
OUT_S2_MOD <- "GSE303006_M23_rev_stage2_energy_by_modality.tsv"
OUT_LOG <- "GSE303006_M23_rev_stage1_stage2_energy_log.txt"

for (f in c(STAGE1_FILES, DIAG_FILE, CORE_FILE, SMALL_FILE)) {
    if (!file.exists(f)) stop("Missing input: ", f)
}

log_con <- rev_start_log(OUT_LOG, "19_stage1_stage2_energy_report.R")

modalities <- names(STAGE1_FILES)

# ------------------------------------------------------------
# 1. Stage 1
# ------------------------------------------------------------

S1 <- lapply(STAGE1_FILES, readRDS)

s1_rows <- list()
cum_rows <- list()

for (m in modalities) {

    s <- S1[[m]]
    d <- as.numeric(s$singular_values)
    ef <- as.numeric(s$explained_fraction)
    K <- length(d)

    # Total energy recovered from 05's definition
    # explained_fraction = d^2 / ||Z||^2
    total_fro2 <- d[1]^2 / ef[1]

    retained_fro2 <- sum(d^2)
    f1 <- retained_fro2 / total_fro2

    s1_rows[[m]] <- data.table(
        modality = m,
        n_features = nrow(s$feature_scores),
        n_cells = nrow(s$cell_loadings),
        K = K,
        total_energy = total_fro2,
        retained_energy = retained_fro2,
        retained_fraction = f1,
        relative_reconstruction_error = sqrt(1 - f1),
        first_component_fraction = ef[1],
        last_component_fraction = ef[K],
        d_K_over_d_1 = d[K] / d[1]
    )

    cum_rows[[m]] <- data.table(
        modality = m,
        k = seq_len(K),
        singular_value = d,
        component_fraction = ef,
        cumulative_fraction = cumsum(ef)
    )
}

s1 <- rbindlist(s1_rows)
s1_cum <- rbindlist(cum_rows)

cat("===== Stage 1: variance retention per modality =====\n")
print(s1[, .(
    modality, n_features, n_cells, K,
    retained_fraction = signif(retained_fraction, 4),
    rel_error = signif(relative_reconstruction_error, 4),
    first_comp = signif(first_component_fraction, 4),
    last_comp = signif(last_component_fraction, 4),
    dK_d1 = signif(d_K_over_d_1, 4)
)])

# ------------------------------------------------------------
# 2. Stage 2
# ------------------------------------------------------------

D <- readRDS(DIAG_FILE)
CORE <- readRDS(CORE_FILE)
SF <- readRDS(SMALL_FILE)

tensor_fro2 <- as.numeric(D$tensor_frobenius2)
core_fro2 <- sum(CORE^2)
modal_norm2 <- D$modality_frobenius2
modal_norm2 <- modal_norm2[modalities]

if (abs(core_fro2 - D$core_frobenius2) > 1e-8 * max(1, core_fro2)) {
    warning("Core norm differs from 08 diagnostics.")
}

f2 <- core_fro2 / tensor_fro2

s2 <- data.table(
    quantity = c(
        "tensor_frobenius2",
        "core_frobenius2",
        "retained_fraction",
        "relative_reconstruction_error",
        "EP_mode_energy_fraction",
        "cell_mode_energy_fraction",
        "modality_mode_energy_fraction"
    ),
    value = c(
        tensor_fro2,
        core_fro2,
        f2,
        sqrt(1 - f2),
        D$EP_mode_energy_fraction,
        D$cell_mode_energy_fraction,
        D$modality_mode_energy_fraction
    )
)

cat("\n===== Stage 2: conceptual tensor =====\n")
cat("Tucker ranks: ",
    paste(names(D$tucker_ranks), D$tucker_ranks, sep = "=", collapse = ", "),
    "\n", sep = "")
print(s2[, .(quantity, value = signif(value, 6))])

# ------------------------------------------------------------
# 3. Stage-2 retention per modality
# ------------------------------------------------------------

U_MODAL <- SF$U_MODAL
U_K <- SF$U_K

stopifnot(
    dim(CORE)[3] == nrow(U_K),
    max(abs(crossprod(U_K) - diag(nrow(U_K)))) < 1e-10,
    dim(CORE)[4] == nrow(U_MODAL),
    max(abs(crossprod(U_MODAL) - diag(ncol(U_MODAL)))) < 1e-10
)

dc <- dim(CORE)
Gpre <- array(0, dim = c(dc[1], dc[2], dc[3], nrow(U_MODAL)))
for (im in seq_len(nrow(U_MODAL))) {
    for (r in seq_len(ncol(U_MODAL))) {
        Gpre[, , , im] <- Gpre[, , , im] + CORE[, , , r] * U_MODAL[im, r]
    }
}

captured <- vapply(seq_along(modalities),
                   function(im) sum(Gpre[, , , im]^2), numeric(1))
names(captured) <- modalities

check_sum <- abs(sum(captured) - core_fro2) / core_fro2
cat("\nCheck: sum of per-modality captured energy vs core norm^2, ",
    "relative difference = ", format(check_sum, scientific = TRUE), "\n",
    sep = "")
if (check_sum > 1e-8) warning("Per-modality decomposition of core energy failed.")

s2_mod <- data.table(
    modality = modalities,
    stage2_slice_energy = as.numeric(modal_norm2),
    mapping_norm_deviation = as.numeric(modal_norm2) - 1,
    stage2_captured_energy = as.numeric(captured),
    stage2_retained_fraction = as.numeric(captured / modal_norm2)
)
s2_mod <- merge(
    s2_mod,
    s1[, .(modality, stage1_retained_fraction = retained_fraction)],
    by = "modality", sort = FALSE
)
s2_mod[, overall_retained_fraction :=
           stage1_retained_fraction * stage2_retained_fraction]

cat("\n===== Retention per modality: Stage 1, Stage 2, overall =====\n")
cat("(stage2_slice_energy should be 1 if the degree-corrected mapping\n",
    " preserves the Frobenius norm of the Stage-1 reconstruction)\n", sep = "")
print(s2_mod[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)])

# ------------------------------------------------------------
# 4. Save
# ------------------------------------------------------------

fwrite(s1, OUT_S1, sep = "\t")
fwrite(s1_cum, OUT_S1_CUM, sep = "\t")
fwrite(s2, OUT_S2, sep = "\t")
fwrite(s2_mod, OUT_S2_MOD, sep = "\t")

cat("\nSaved:\n", OUT_S1, "\n", OUT_S1_CUM, "\n", OUT_S2, "\n",
    OUT_S2_MOD, "\n", OUT_LOG, "\n", sep = "")

rev_end_log(log_con)
