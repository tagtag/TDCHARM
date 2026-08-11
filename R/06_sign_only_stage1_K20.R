# ============================================================
# 06_sign_only_stage1_K20.R
# GSE303006 / GENCODE M23
#
# Sign-only orientation of independently obtained K=20 SVD components.
#
# IMPORTANT:
# - No Procrustes rotation.
# - No component mixing.
# - No permutation of components.
# - No V_common.
#
# For each modality and each component k:
#   find the cell with the largest absolute cell loading.
#   If that loading is negative, multiply BOTH
#       feature_scores[, k]
#       cell_loadings[, k]
#   by -1.
#
# This fixes a deterministic orientation while preserving:
#
#   feature_scores[,k] %*% t(cell_loadings[,k])
#
# exactly.
#
# Stage-2 tensor remains conceptually:
#   730,969 E-P x 4,258 cells x 20 components x 4 modalities
# ============================================================

library(data.table)

K <- 20L

INPUT_FILES <- c(
    ATAC     = "GSE303006_M23_stage1_K20_ATAC.rds",
    H3K27ac  = "GSE303006_M23_stage1_K20_H3K27ac.rds",
    RNA      = "GSE303006_M23_stage1_K20_RNA.rds",
    `3D`     = "GSE303006_M23_stage1_K20_3D.rds"
)

OUTPUT_FILES <- c(
    ATAC     = "GSE303006_M23_stage1_K20_ATAC_signaligned.rds",
    H3K27ac  = "GSE303006_M23_stage1_K20_H3K27ac_signaligned.rds",
    RNA      = "GSE303006_M23_stage1_K20_RNA_signaligned.rds",
    `3D`     = "GSE303006_M23_stage1_K20_3D_signaligned.rds"
)

DIAG_FILE <- "GSE303006_M23_stage1_K20_sign_alignment_diagnostics.tsv"

missing <- INPUT_FILES[!file.exists(INPUT_FILES)]
if (length(missing) > 0L) {
    stop(
        "Missing Stage-1 file(s):\n",
        paste(missing, collapse = "\n")
    )
}

sign_align_one <- function(obj, modality, K = 20L) {

    required <- c(
        "feature_scores",
        "cell_loadings",
        "singular_values"
    )

    if (!all(required %in% names(obj))) {
        stop(
            modality,
            ": missing object field(s): ",
            paste(setdiff(required, names(obj)), collapse = ", ")
        )
    }

    B <- obj$feature_scores
    V <- obj$cell_loadings

    if (ncol(B) != K || ncol(V) != K) {
        stop(
            modality,
            ": expected K=", K,
            ", but dimensions are B=",
            paste(dim(B), collapse = "x"),
            ", V=",
            paste(dim(V), collapse = "x")
        )
    }

    signs <- rep(1, K)
    anchor_cell <- character(K)
    anchor_before <- numeric(K)

    # A small subset of cells is sufficient for numerical
    # reconstruction-preservation checking.
    check_rows <- seq_len(min(25L, nrow(V)))

    B_before <- B
    V_before <- V

    for (k in seq_len(K)) {

        j <- which.max(abs(V[, k]))

        anchor_cell[k] <- rownames(V)[j]
        anchor_before[k] <- V[j, k]

        if (V[j, k] < 0) {

            signs[k] <- -1

            V[, k] <- -V[, k]
            B[, k] <- -B[, k]
        }
    }

    # Check that each rank-1 reconstruction is unchanged.
    preservation_error <- numeric(K)

    for (k in seq_len(K)) {

        old_block <- tcrossprod(
            B_before[, k],
            V_before[check_rows, k]
        )

        new_block <- tcrossprod(
            B[, k],
            V[check_rows, k]
        )

        preservation_error[k] <- max(
            abs(old_block - new_block)
        )
    }

    # Confirm every anchor is now positive.
    anchor_after <- numeric(K)

    for (k in seq_len(K)) {

        j <- which.max(abs(V[, k]))
        anchor_after[k] <- V[j, k]

        if (anchor_after[k] < 0) {
            stop(
                modality,
                " component ", k,
                ": sign convention failed."
            )
        }
    }

    obj$feature_scores <- B
    obj$cell_loadings <- V
    obj$sign_alignment <- list(
        method = "largest-absolute-cell-loading-positive",
        signs = signs,
        anchor_cell = anchor_cell,
        anchor_before = anchor_before,
        anchor_after = anchor_after,
        reconstruction_preservation_error = preservation_error
    )

    diag <- data.table(
        modality = modality,
        component = seq_len(K),
        sign_applied = signs,
        flipped = signs == -1,
        anchor_cell = anchor_cell,
        anchor_loading_before = anchor_before,
        anchor_loading_after = anchor_after,
        reconstruction_preservation_error = preservation_error
    )

    list(
        object = obj,
        diagnostic = diag
    )
}


all_diag <- list()

for (nm in names(INPUT_FILES)) {

    cat("\nProcessing ", nm, "...\n", sep = "")

    obj <- readRDS(INPUT_FILES[[nm]])

    ans <- sign_align_one(
        obj = obj,
        modality = nm,
        K = K
    )

    saveRDS(
        ans$object,
        OUTPUT_FILES[[nm]]
    )

    all_diag[[nm]] <- ans$diagnostic

    cat(
        "  flipped components : ",
        sum(ans$diagnostic$flipped),
        " / ", K, "\n",
        sep = ""
    )

    cat(
        "  max reconstruction preservation error : ",
        format(
            max(ans$diagnostic$reconstruction_preservation_error),
            scientific = TRUE
        ),
        "\n",
        sep = ""
    )
}

diag_dt <- rbindlist(
    all_diag,
    use.names = TRUE
)

fwrite(
    diag_dt,
    DIAG_FILE,
    sep = "\t"
)

cat("\n============================================\n")
cat("Sign-only alignment completed\n")
cat("============================================\n")

cat("\nOutput files:\n")
cat(paste(OUTPUT_FILES, collapse = "\n"), "\n")
cat(DIAG_FILE, "\n")

cat("\nSign flips by modality:\n")
print(
    diag_dt[
        ,
        .(
            flipped = sum(flipped),
            unchanged = sum(!flipped)
        ),
        by = modality
    ]
)

cat(
    "\nStage-2 conceptual tensor remains:\n",
    "730969 x 4258 x 20 x 4\n",
    sep = ""
)

cat(
    "\nNote: sign-only alignment changes interpretation/orientation only.\n",
    "Because feature_scores and cell_loadings are flipped together,\n",
    "their rank-1 product is unchanged.\n",
    sep = ""
)
