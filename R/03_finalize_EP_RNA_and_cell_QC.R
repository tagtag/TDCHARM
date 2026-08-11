# ============================================================
# 03_finalize_EP_RNA_and_cell_QC.R
# GSE303006 / GENCODE M23
#
# This script fills the reproducibility gap between:
#   02_build_M23_EP_and_3D_summary.R
# and
#   04_build_3D_HDF5.R / 05_stage1_independent_K20.R
#
# It performs the final filters used in the manuscript:
#
# 1. 3D bin-pair QC:
#      coverage_3D >= 0.90
#      finite, non-zero 3D distance SD
#
# 2. Enhancer QC across the 4,383 common ATAC/H3K27ac/3D cells:
#      ATAC detection fraction >= 0.005
#      H3K27ac detection fraction >= 0.005
#
# 3. RNA extraction for the 16,239 final target genes:
#      library size calculated from ALL RNA genes
#      selected-gene matrix normalized as logCP10K
#
# 4. Cell QC:
#      ATAC enhancer-associated fragments >= 250
#      H3K27ac enhancer-associated fragments >= 100
#      RNA library size >= 500
#      RNA detected selected genes >= 500
#
# Expected reported final dimensions:
#      730,969 E-P pairs
#       42,669 enhancers
#       16,239 genes
#      391,435 20-kb bin pairs
#        4,258 final cells
#
# Linux/Unix requirements:
#      zcat and awk
# ============================================================

library(data.table)
library(Matrix)

# ------------------------------------------------------------
# 0. Input / output
# ------------------------------------------------------------

EP_WITH_BP_FILE <-
    "GSE303006_M23_brain_EP_with_20kb_binpair.rds"

BP_3D_SUMMARY_FILE <-
    "GSE303006_M23_brain_20kb_3D_summary_complete.rds"

ATAC_FILE <-
    "GSE303006_M23_brain_ATAC_enhancer_by_cell_AH3D.rds"

H3_FILE <-
    "GSE303006_M23_brain_H3K27ac_enhancer_by_cell_AH3D.rds"

COMMON_CELLS_FILE <-
    "GSE303006_M23_brain_common_cells_ATAC_H3K27ac_3D.rds"

RNA_COUNTS_FILE <-
    "GSE303006_brain.RNAcounts.gene.total.format.tsv.gz"

OUT_EP <-
    "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"

OUT_EP_TSV <-
    "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.tsv.gz"

OUT_ENH_QC <-
    "GSE303006_M23_brain_enhancer_QC_detection005.tsv.gz"

OUT_RNA_RAW <-
    "GSE303006_M23_brain_RNA_raw_gene_by_cell.rds"

OUT_RNA_LOG <-
    "GSE303006_M23_brain_RNA_logCP10K_gene_by_cell.rds"

OUT_RNA_QC <-
    "GSE303006_M23_brain_RNA_cell_QC.tsv"

OUT_FINAL_QC <-
    "GSE303006_M23_brain_final_cell_QC.tsv"

OUT_FINAL_CELLS <-
    "GSE303006_M23_brain_final_4258_cells.rds"

# Final thresholds used in the manuscript.
MIN_3D_COVERAGE <- 0.90
MIN_ENH_DETECT_FRAC <- 0.005

MIN_ATAC_CELL_COUNT <- 250
MIN_H3_CELL_COUNT <- 100
MIN_RNA_LIBRARY_SIZE <- 500
MIN_RNA_DETECTED_GENES <- 500

SCALE_FACTOR <- 10000

# Save the raw 16,239 x 4,383 RNA sparse matrix.
# Set FALSE if disk space is limited; Stage 1 needs only OUT_RNA_LOG.
SAVE_RNA_RAW <- TRUE

# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

need <- c(
    EP_WITH_BP_FILE,
    BP_3D_SUMMARY_FILE,
    ATAC_FILE,
    H3_FILE,
    COMMON_CELLS_FILE,
    RNA_COUNTS_FILE
)

missing <- need[!file.exists(need)]

if (length(missing) > 0L) {
    stop(
        "Missing input file(s):\n",
        paste(missing, collapse = "\n")
    )
}

if (Sys.which("zcat") == "") {
    stop("zcat is required by this preprocessing script.")
}

if (Sys.which("awk") == "") {
    stop("awk is required by this preprocessing script.")
}

# ------------------------------------------------------------
# 2. Load E-P / 3D summary / enhancer matrices
# ------------------------------------------------------------

ep <- as.data.table(
    readRDS(EP_WITH_BP_FILE)
)

bp3d <- as.data.table(
    readRDS(BP_3D_SUMMARY_FILE)
)

ATAC <- readRDS(ATAC_FILE)
H3 <- readRDS(H3_FILE)

common_cells <- as.character(
    readRDS(COMMON_CELLS_FILE)
)

if (!all(common_cells %in% colnames(ATAC)) ||
    !all(common_cells %in% colnames(H3))) {
    stop("Common-cell list is inconsistent with ATAC/H3K27ac matrices.")
}

ATAC <- ATAC[
    ,
    common_cells,
    drop = FALSE
]

H3 <- H3[
    ,
    common_cells,
    drop = FALSE
]

if (!identical(colnames(ATAC), colnames(H3))) {
    stop("ATAC and H3K27ac cell order differs.")
}

cat(
    "Common ATAC/H3K27ac/3D cells :",
    length(common_cells),
    "\n"
)

# ------------------------------------------------------------
# 3. 3D bin-pair QC
# ------------------------------------------------------------

required_bp <- c(
    "binpair_id",
    "coverage_3D",
    "sd_distance_3D"
)

if (!all(required_bp %in% names(bp3d))) {
    stop(
        "3D summary lacks: ",
        paste(
            setdiff(required_bp, names(bp3d)),
            collapse = ", "
        )
    )
}

usable_bp <- bp3d[
    coverage_3D >= MIN_3D_COVERAGE &
    is.finite(sd_distance_3D) &
    sd_distance_3D > 0,
    as.integer(binpair_id)
]

cat(
    "Usable 3D bin pairs :",
    length(usable_bp),
    "/",
    nrow(bp3d),
    "\n"
)

if (!"same_20k_bin" %in% names(ep)) {
    stop("E-P table lacks same_20k_bin.")
}

ep_3d <- ep[
    same_20k_bin == FALSE &
    binpair_id %in% usable_bp
]

cat(
    "E-P pairs after 3D QC :",
    nrow(ep_3d),
    "\n"
)

# ------------------------------------------------------------
# 4. Enhancer detection QC (0.5% in both ATAC and H3K27ac)
# ------------------------------------------------------------

if (is.null(rownames(ATAC)) ||
    is.null(rownames(H3))) {
    stop("ATAC/H3K27ac matrices require enhancer row names.")
}

if (!identical(rownames(ATAC), rownames(H3))) {
    stop("ATAC/H3K27ac enhancer row order differs.")
}

atac_detect_frac <- as.numeric(
    Matrix::rowMeans(
        ATAC > 0
    )
)

h3_detect_frac <- as.numeric(
    Matrix::rowMeans(
        H3 > 0
    )
)

enh_qc <- data.table(
    enhancer_id = rownames(ATAC),
    ATAC_detection_fraction = atac_detect_frac,
    H3K27ac_detection_fraction = h3_detect_frac
)

enh_qc[
    ,
    keep_enhancer :=
        ATAC_detection_fraction >= MIN_ENH_DETECT_FRAC &
        H3K27ac_detection_fraction >= MIN_ENH_DETECT_FRAC
]

fwrite(
    enh_qc,
    OUT_ENH_QC,
    sep = "\t"
)

keep_enhancers <- enh_qc[
    keep_enhancer == TRUE,
    enhancer_id
]

ep_final <- ep_3d[
    enhancer_id %in% keep_enhancers
]

# Ensure the final E-P pair definition is unique.
dup <- ep_final[
    ,
    .N,
    by = .(
        enhancer_id,
        gene_name
    )
][
    N > 1L
]

if (nrow(dup) > 0L) {
    stop(
        "Duplicate enhancer_id x gene_name pairs remain: ",
        nrow(dup)
    )
}

saveRDS(
    ep_final,
    OUT_EP
)

fwrite(
    ep_final,
    OUT_EP_TSV,
    sep = "\t"
)

genes <- sort(
    unique(
        as.character(
            ep_final$gene_name
        )
    )
)

enhancers <- sort(
    unique(
        as.character(
            ep_final$enhancer_id
        )
    )
)

binpairs <- sort(
    unique(
        as.integer(
            ep_final$binpair_id
        )
    )
)

cat(
    "\nFinal E-P pairs :",
    nrow(ep_final),
    "\n"
)

cat(
    "Unique enhancers :",
    length(enhancers),
    "\n"
)

cat(
    "Unique genes :",
    length(genes),
    "\n"
)

cat(
    "Unique bin-pairs :",
    length(binpairs),
    "\n"
)

# ------------------------------------------------------------
# 5. RNA library sizes using ALL RNA genes
# ------------------------------------------------------------

cat(
    "\nCalculating RNA library sizes using all genes...\n"
)

awk_lib <- paste0(
    "NR==1 {",
    "for(i=2;i<=NF;i++){name[i]=$i}; n=NF; next",
    "} ",
    "{for(i=2;i<=NF;i++){s[i]+=$i}} ",
    "END {for(i=2;i<=n;i++){print name[i] \"\\t\" s[i]}}"
)

cmd_lib <- sprintf(
    "zcat %s | awk -F '\\t' '%s'",
    shQuote(RNA_COUNTS_FILE),
    awk_lib
)

rna_lib <- fread(
    cmd = cmd_lib,
    col.names = c(
        "cell",
        "library_size"
    )
)

if (!all(common_cells %in% rna_lib$cell)) {
    stop("Some common cells are absent from the RNA count matrix.")
}

lib_size <- rna_lib$library_size[
    match(
        common_cells,
        rna_lib$cell
    )
]

names(lib_size) <- common_cells

cat(
    "Library size summary:\n"
)

print(
    summary(
        lib_size
    )
)

# ------------------------------------------------------------
# 6. Extract only final 16,239 genes from the RNA count matrix
# ------------------------------------------------------------

cat(
    "\nReading selected RNA genes...\n"
)

gene_file <- tempfile(
    fileext = ".txt"
)

writeLines(
    genes,
    gene_file
)

# First awk input is gene_file; second input is stdin ("-").
# The RNA header is retained, and only requested gene rows are emitted.
awk_select <- paste0(
    "NR==FNR {keep[$1]=1; next} ",
    "FNR==1 || ($1 in keep)"
)

cmd_select <- sprintf(
    "zcat %s | awk -F '\\t' '%s' %s -",
    shQuote(RNA_COUNTS_FILE),
    awk_select,
    shQuote(gene_file)
)

rna_dt <- fread(
    cmd = cmd_select,
    check.names = FALSE
)

unlink(
    gene_file
)

gene_col <- names(rna_dt)[1]

setnames(
    rna_dt,
    gene_col,
    "gene_name"
)

missing_genes <- setdiff(
    genes,
    rna_dt$gene_name
)

cat(
    "RNA rows extracted :",
    nrow(rna_dt),
    "\n"
)

cat(
    "Requested genes absent from RNA matrix :",
    length(missing_genes),
    "\n"
)

if (length(missing_genes) > 0L) {
    stop(
        "Missing requested RNA genes: ",
        paste(
            head(missing_genes, 20),
            collapse = ", "
        )
    )
}

dup_gene <- rna_dt[
    ,
    .N,
    by = gene_name
][
    N > 1L
]

cat(
    "Duplicated gene symbols in RNA matrix :",
    nrow(dup_gene),
    "\n"
)

if (nrow(dup_gene) > 0L) {
    stop(
        "Duplicated requested gene symbols found."
    )
}

if (!all(common_cells %in% names(rna_dt))) {
    stop(
        "Some common cells are absent from extracted RNA columns."
    )
}

rna_dt <- rna_dt[
    match(
        genes,
        gene_name
    )
]

rna_dense <- as.matrix(
    rna_dt[
        ,
        ..common_cells
    ]
)

storage.mode(
    rna_dense
) <- "numeric"

rna_raw <- Matrix::Matrix(
    rna_dense,
    sparse = TRUE
)

rna_raw <- as(
    rna_raw,
    "dgCMatrix"
)

rownames(
    rna_raw
) <- genes

colnames(
    rna_raw
) <- common_cells

rm(
    rna_dense,
    rna_dt
)

gc()

cat(
    "RNA raw dimensions :",
    paste(
        dim(rna_raw),
        collapse = " x "
    ),
    "\n"
)

cat(
    "RNA raw nnzero :",
    Matrix::nnzero(rna_raw),
    "\n"
)

if (SAVE_RNA_RAW) {
    saveRDS(
        rna_raw,
        OUT_RNA_RAW
    )
}

# ------------------------------------------------------------
# 7. RNA QC and logCP10K
# ------------------------------------------------------------

detected_genes <- as.numeric(
    Matrix::colSums(
        rna_raw > 0
    )
)

names(
    detected_genes
) <- common_cells

rna_qc <- data.table(
    cell = common_cells,
    library_size = as.numeric(
        lib_size[
            common_cells
        ]
    ),
    detected_genes = as.numeric(
        detected_genes[
            common_cells
        ]
    )
)

rna_qc[
    ,
    keep_RNA :=
        library_size >= MIN_RNA_LIBRARY_SIZE &
        detected_genes >= MIN_RNA_DETECTED_GENES
]

fwrite(
    rna_qc,
    OUT_RNA_QC,
    sep = "\t"
)

cat(
    "\nRNA QC pass :",
    sum(rna_qc$keep_RNA),
    "/",
    nrow(rna_qc),
    "\n"
)

if (any(lib_size <= 0)) {
    stop("RNA library size <= 0 encountered.")
}

rna_log <- rna_raw %*%
    Matrix::Diagonal(
        x =
            SCALE_FACTOR /
            as.numeric(
                lib_size[
                    common_cells
                ]
            )
    )

rna_log <- as(
    rna_log,
    "dgCMatrix"
)

rna_log@x <- log1p(
    rna_log@x
)

rownames(
    rna_log
) <- genes

colnames(
    rna_log
) <- common_cells

saveRDS(
    rna_log,
    OUT_RNA_LOG
)

# ------------------------------------------------------------
# 8. Final ATAC / H3K27ac / RNA cell QC
# ------------------------------------------------------------

atac_cell <- as.numeric(
    Matrix::colSums(
        ATAC
    )
)

h3_cell <- as.numeric(
    Matrix::colSums(
        H3
    )
)

names(atac_cell) <- common_cells
names(h3_cell) <- common_cells

final_qc <- data.table(
    cell = common_cells,
    ATAC_count = atac_cell,
    H3K27ac_count = h3_cell,
    RNA_library_size = as.numeric(
        lib_size[
            common_cells
        ]
    ),
    RNA_detected_genes = as.numeric(
        detected_genes[
            common_cells
        ]
    )
)

final_qc[
    ,
    keep_ATAC :=
        ATAC_count >=
        MIN_ATAC_CELL_COUNT
]

final_qc[
    ,
    keep_H3K27ac :=
        H3K27ac_count >=
        MIN_H3_CELL_COUNT
]

final_qc[
    ,
    keep_RNA :=
        RNA_library_size >=
        MIN_RNA_LIBRARY_SIZE &
        RNA_detected_genes >=
        MIN_RNA_DETECTED_GENES
]

# All rows already belong to the common ATAC/H3K27ac/3D set.
final_qc[
    ,
    keep_3D := TRUE
]

final_qc[
    ,
    keep_final :=
        keep_ATAC &
        keep_H3K27ac &
        keep_RNA &
        keep_3D
]

final_cells <- final_qc[
    keep_final == TRUE,
    cell
]

fwrite(
    final_qc,
    OUT_FINAL_QC,
    sep = "\t"
)

saveRDS(
    final_cells,
    OUT_FINAL_CELLS
)

cat(
    "\nFinal 4-modality cells :",
    length(final_cells),
    "\n"
)

# ------------------------------------------------------------
# 9. Reproducibility checks against reported article dimensions
# ------------------------------------------------------------

expected <- list(
    EP_pairs = 730969L,
    enhancers = 42669L,
    genes = 16239L,
    binpairs = 391435L,
    cells = 4258L
)

observed <- list(
    EP_pairs = nrow(ep_final),
    enhancers = length(enhancers),
    genes = length(genes),
    binpairs = length(binpairs),
    cells = length(final_cells)
)

cat(
    "\nReported-dimension checks:\n"
)

for (nm in names(expected)) {

    ok <- identical(
        as.integer(observed[[nm]]),
        as.integer(expected[[nm]])
    )

    cat(
        sprintf(
            "%-10s observed=%d expected=%d %s\n",
            nm,
            observed[[nm]],
            expected[[nm]],
            if (ok) "OK" else "WARNING"
        )
    )

    if (!ok) {
        warning(
            nm,
            ": observed ",
            observed[[nm]],
            " differs from reported ",
            expected[[nm]],
            "."
        )
    }
}

cat(
    "\n03_finalize_EP_RNA_and_cell_QC.R finished.\n"
)
