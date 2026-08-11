# ============================================================
# 03_build_3D_HDF5.R
# GSE303006 / GENCODE M23
# Final 20-kb bin-pair x cell 3D-distance matrix -> HDF5
# Resume-capable, float32, block-wise
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE))
    stop("Install data.table first.")
if (!requireNamespace("rhdf5", quietly = TRUE))
    stop("Install rhdf5 with BiocManager::install('rhdf5').")

library(data.table)
library(rhdf5)

# ---------- settings ----------
ep_file    <- "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"
bp_file    <- "GSE303006_M23_brain_20kb_binpairs_indexed.rds"
loci_file  <- "GSE303006_M23_brain_20kb_loci_for_EP.rds"
cells_file <- "GSE303006_M23_brain_common_cells_ATAC_H3K27ac_3D.rds"

tdg_dir <- "tdg20k/geo_submit/brain_20k_tdg"

bp_final_file <- "GSE303006_M23_brain_binpairs_final_enhQC005.rds"
h5file <- "GSE303006_M23_brain_3D_binpair_by_cell.h5"

BLOCK_CELLS <- 16L

# TRUE にすると既存HDF5を削除して最初から作る。
# 通常は FALSE。
RESET_H5 <- FALSE

# ---------- preflight ----------
required_files <- c(ep_file, bp_file, loci_file, cells_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L)
    stop("Missing files:\n", paste(missing_files, collapse = "\n"))

if (!dir.exists(tdg_dir))
    stop("3D directory not found: ", tdg_dir)

if (RESET_H5 && file.exists(h5file)) {
    h5closeAll()
    unlink(h5file)
    cat("Existing HDF5 removed.\n")
}

# ---------- load ----------
ep_final_qc <- as.data.table(readRDS(ep_file))
bp_all      <- as.data.table(readRDS(bp_file))
loci        <- as.data.table(readRDS(loci_file))
cells_run   <- as.character(readRDS(cells_file))

need_bp <- c("binpair_id", "enhancer_chr", "enhancer_bin",
             "promoter_bin", "e_locus", "p_locus")
need_loci <- c("locus_id", "key_mat", "key_pat")

if (!all("binpair_id" %in% names(ep_final_qc)))
    stop("binpair_id missing from ep_final_qc.")
if (!all(need_bp %in% names(bp_all)))
    stop("Required columns missing from bp_all.")
if (!all(need_loci %in% names(loci)))
    stop("Required columns missing from loci.")

# ---------- final unique bin-pairs ----------
binpair_keep <- sort(unique(as.integer(ep_final_qc$binpair_id)))

idx <- match(binpair_keep, as.integer(bp_all$binpair_id))
if (anyNA(idx))
    stop("Some final binpair_id values are absent from bp_all.")

bp_final <- bp_all[idx]

if (!identical(as.integer(bp_final$binpair_id), binpair_keep))
    stop("binpair order mismatch.")

if (anyNA(bp_final$e_locus) || anyNA(bp_final$p_locus))
    stop("NA in e_locus/p_locus.")

if (any(bp_final$e_locus < 1L) ||
    any(bp_final$p_locus < 1L) ||
    any(bp_final$e_locus > nrow(loci)) ||
    any(bp_final$p_locus > nrow(loci)))
    stop("Out-of-range locus index.")

saveRDS(bp_final, bp_final_file)

cat("Final E-P pairs :", nrow(ep_final_qc), "\n")
cat("Final bin-pairs :", nrow(bp_final), "\n")
cat("Cells           :", length(cells_run), "\n")

# ---------- 3dg files ----------
tdg_files <- list.files(
    tdg_dir,
    pattern = "\\.20k\\.3dg\\.gz$",
    full.names = TRUE
)

tdg_cells <- sub("\\.20k\\.3dg\\.gz$", "", basename(tdg_files))

if (anyDuplicated(tdg_cells))
    stop("Duplicated 3dg cell names.")

file_by_cell <- setNames(tdg_files, tdg_cells)

missing_3d <- setdiff(cells_run, names(file_by_cell))
if (length(missing_3d) > 0L)
    stop("Missing 3dg cells:\n", paste(missing_3d, collapse = " "))

# ---------- read 3dg ----------
read_3dg <- function(file) {

    d <- fread(
        cmd = paste("zcat", shQuote(file)),
        header = FALSE,
        col.names = c("chrom_phase", "position", "x", "y", "z")
    )

    d[, chr := sub("\\((mat|pat)\\)$", "", chrom_phase)]
    d[, phase := sub("^.*\\((mat|pat)\\)$", "\\1", chrom_phase)]
    d[, key := paste(chr, position, phase, sep = ":")]

    d
}

# ---------- one-cell distance ----------
# cis E-P: mat-mat and pat-pat only.
# If both exist, average them. If only one exists, use the available one.
calc_3D_fast <- function(bp, loci, tdg) {

    idx_mat <- match(loci$key_mat, tdg$key)
    idx_pat <- match(loci$key_pat, tdg$key)

    em <- idx_mat[bp$e_locus]
    pm <- idx_mat[bp$p_locus]

    epat <- idx_pat[bp$e_locus]
    ppat <- idx_pat[bp$p_locus]

    n <- nrow(bp)

    d_mat <- rep(NA_real_, n)
    d_pat <- rep(NA_real_, n)

    ok_mat <- !is.na(em) & !is.na(pm)

    if (any(ok_mat)) {
        d_mat[ok_mat] <- sqrt(
            (tdg$x[em[ok_mat]] - tdg$x[pm[ok_mat]])^2 +
            (tdg$y[em[ok_mat]] - tdg$y[pm[ok_mat]])^2 +
            (tdg$z[em[ok_mat]] - tdg$z[pm[ok_mat]])^2
        )
    }

    ok_pat <- !is.na(epat) & !is.na(ppat)

    if (any(ok_pat)) {
        d_pat[ok_pat] <- sqrt(
            (tdg$x[epat[ok_pat]] - tdg$x[ppat[ok_pat]])^2 +
            (tdg$y[epat[ok_pat]] - tdg$y[ppat[ok_pat]])^2 +
            (tdg$z[epat[ok_pat]] - tdg$z[ppat[ok_pat]])^2
        )
    }

    d_mean <- rep(NA_real_, n)

    both <- !is.na(d_mat) & !is.na(d_pat)
    only_mat <- !is.na(d_mat) & is.na(d_pat)
    only_pat <- is.na(d_mat) & !is.na(d_pat)

    d_mean[both] <- (d_mat[both] + d_pat[both]) / 2
    d_mean[only_mat] <- d_mat[only_mat]
    d_mean[only_pat] <- d_pat[only_pat]

    d_mean
}

# ---------- test one cell ----------
test_cell <- cells_run[1]

cat("\nTesting one cell:", test_cell, "\n")

tdg_test <- read_3dg(file_by_cell[[test_cell]])
v_test <- calc_3D_fast(bp_final, loci, tdg_test)

test_cov <- mean(!is.na(v_test))

cat("Test coverage :", test_cov, "\n")
print(summary(v_test))

if (!is.finite(test_cov) || test_cov < 0.5)
    stop("Unexpectedly low 3D coverage; stopping before HDF5 writing.")

rm(tdg_test, v_test)
gc()

# ---------- HDF5 dimensions ----------
n_bp   <- nrow(bp_final)
n_cell <- length(cells_run)

float32_type <- "H5T_IEEE_F32LE"

# h5createDataset() accepts an explicit HDF5 type name.
if (!(float32_type %in% rhdf5::h5const("H5T")))
    stop("H5T_IEEE_F32LE is not available in this rhdf5 build.")

chunk_rows <- min(8192L, n_bp)
chunk_cols <- min(BLOCK_CELLS, n_cell)

cat("\nHDF5 dimensions :", n_bp, "x", n_cell, "\n")
cat("HDF5 datatype   : float32\n")
cat("HDF5 chunk      :", chunk_rows, "x", chunk_cols, "\n")

# ---------- create or resume HDF5 ----------
if (!file.exists(h5file)) {

    cat("\nCreating:", h5file, "\n")

    h5createFile(h5file)

    h5createDataset(
        file = h5file,
        dataset = "distance",
        dims = c(n_bp, n_cell),
        H5type = float32_type,
        chunk = c(chunk_rows, chunk_cols),
        level = 4,
        filter = "gzip",
        shuffle = TRUE
    )

    h5createDataset(
        file = h5file,
        dataset = "processed",
        dims = n_cell,
        storage.mode = "integer",
        chunk = min(1024L, n_cell),
        level = 0
    )

    h5write(
        integer(n_cell),
        file = h5file,
        name = "processed"
    )

    h5write(
        as.integer(bp_final$binpair_id),
        file = h5file,
        name = "binpair_id"
    )

    h5write(
        cells_run,
        file = h5file,
        name = "cell"
    )

    h5write(
        as.integer(c(n_bp, n_cell)),
        file = h5file,
        name = "matrix_dim"
    )

    h5write(
        "mean of available cis mat-mat and pat-pat 3D distances",
        file = h5file,
        name = "distance_definition"
    )

    h5closeAll()

    cat("HDF5 initialized.\n")

} else {

    cat("\nExisting HDF5 found: resume mode\n")

    contents <- h5ls(h5file)
    root_names <- contents$name[contents$group == "/"]

    required_ds <- c(
        "distance",
        "processed",
        "binpair_id",
        "cell",
        "matrix_dim"
    )

    missing_ds <- setdiff(required_ds, root_names)

    if (length(missing_ds) > 0L) {
        stop(
            "Existing HDF5 is incomplete. Missing: ",
            paste(missing_ds, collapse = ", "),
            "\nSet RESET_H5 <- TRUE and rerun if this is a failed initialization."
        )
    }
}

# ---------- metadata validation ----------
stored_bp <- as.integer(h5read(h5file, "binpair_id"))
stored_cells <- as.character(h5read(h5file, "cell"))
stored_dim <- as.integer(h5read(h5file, "matrix_dim"))

if (!identical(stored_bp, as.integer(bp_final$binpair_id)))
    stop("Stored binpair_id differs from current bp_final.")

if (!identical(stored_cells, cells_run))
    stop("Stored cell order differs from current cells_run.")

if (!identical(stored_dim, as.integer(c(n_bp, n_cell))))
    stop("Stored matrix dimensions differ from current dimensions.")

processed <- as.integer(h5read(h5file, "processed"))

if (length(processed) != n_cell ||
    any(!processed %in% c(0L, 1L)))
    stop("Invalid processed vector.")

cat("\nAlready processed :", sum(processed == 1L), "/", n_cell, "\n")
cat("Remaining         :", sum(processed == 0L), "\n")

# ---------- main block-wise processing ----------
todo <- which(processed == 0L)

if (length(todo) == 0L) {

    cat("\nAll cells are already processed.\n")

} else {

    blocks <- split(
        todo,
        ceiling(seq_along(todo) / BLOCK_CELLS)
    )

    cat("\nStarting/resuming HDF5 writing\n")
    cat("Remaining cells :", length(todo), "\n")
    cat("Blocks          :", length(blocks), "\n")

    for (b in seq_along(blocks)) {

        cols <- as.integer(blocks[[b]])

        block <- matrix(
            NA_real_,
            nrow = n_bp,
            ncol = length(cols)
        )

        for (j in seq_along(cols)) {

            k <- cols[j]
            cell <- cells_run[k]

            tdg <- read_3dg(
                file_by_cell[[cell]]
            )

            block[, j] <- calc_3D_fast(
                bp_final,
                loci,
                tdg
            )

            rm(tdg)
        }

        # Write distances first.
        h5write(
            block,
            file = h5file,
            name = "distance",
            index = list(NULL, cols)
        )

        # Mark complete only after distance write succeeds.
        h5write(
            rep(1L, length(cols)),
            file = h5file,
            name = "processed",
            index = list(cols)
        )

        h5closeAll()

        cov <- colMeans(!is.na(block))

        cat(
            "block ", b, "/", length(blocks),
            " | cells ", min(cols), "-", max(cols),
            " | median coverage = ",
            sprintf("%.4f", median(cov)),
            "\n",
            sep = ""
        )

        rm(block)

        if (b %% 10L == 0L)
            gc()
    }
}

# ---------- completion check ----------
processed <- as.integer(h5read(h5file, "processed"))

n_done <- sum(processed == 1L)
n_remaining <- sum(processed == 0L)

cat("\nProcessed :", n_done, "/", length(processed), "\n")
cat("Remaining :", n_remaining, "\n")

if (n_remaining == 0L) {
    cat("3D HDF5 matrix is complete.\n")
} else {
    cat("Incomplete but resumable: rerun this same script.\n")
}

# ---------- small read-back validation ----------
test <- h5read(
    h5file,
    "distance",
    index = list(
        seq_len(min(10L, n_bp)),
        seq_len(min(5L, n_cell))
    )
)

cat("\nSmall read-back test:\n")
print(test)

cat("\nRead-back summary:\n")
print(summary(as.numeric(test)))

cat("\nHDF5 contents:\n")
print(h5ls(h5file))

h5closeAll()

# ---------- file size ----------
size_gib <- file.info(h5file)$size / 1024^3

cat("\nHDF5 file size (GiB):", round(size_gib, 3), "\n")
cat("\n03_build_3D_HDF5.R finished.\n")
