# ============================================================
# 07_prepare_stage2_implicit.R
# GSE303006 / GENCODE M23
#
# Prepare the Stage-2 hierarchical Tucker input WITHOUT
# materializing the huge 4-way tensor.
#
# Conceptual tensor:
#
#   X[p, c, k, m]
#   730,969 E-P pairs x 4,258 cells x 20 components x 4 modalities
#
# with
#
#   ATAC:
#     X[p,c,k,ATAC] =
#       B_ATAC[e(p),k] * V_ATAC[c,k] / sqrt(d_e)
#
#   H3K27ac:
#     X[p,c,k,H3K27ac] =
#       B_H3[e(p),k] * V_H3[c,k] / sqrt(d_e)
#
#   RNA:
#     X[p,c,k,RNA] =
#       B_RNA[g(p),k] * V_RNA[c,k] / sqrt(d_g)
#
#   3D:
#     X[p,c,k,3D] =
#       B_3D[b(p),k] * V_3D[c,k] / sqrt(d_b)
#
# where:
#   d_e = number of final E-P pairs containing enhancer e
#   d_g = number of final E-P pairs containing gene g
#   d_b = number of final E-P pairs sharing 20-kb bin-pair b
#
# The 1/sqrt(d) correction guarantees that replication of one
# underlying enhancer/gene/bin-pair over multiple E-P rows does
# not inflate its total squared contribution.
#
# This script creates:
#   1) E-P mapping + degree table
#   2) a compact manifest for the implicit tensor
#   3) helper functions to generate only requested tensor blocks
#
# It does NOT run Stage-2 Tucker decomposition yet.
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
}

library(data.table)

K <- 20L

EP_FILE <-
    "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"

STAGE1_FILES <- c(
    ATAC =
        "GSE303006_M23_stage1_K20_ATAC_signaligned.rds",
    H3K27ac =
        "GSE303006_M23_stage1_K20_H3K27ac_signaligned.rds",
    RNA =
        "GSE303006_M23_stage1_K20_RNA_signaligned.rds",
    `3D` =
        "GSE303006_M23_stage1_K20_3D_signaligned.rds"
)

MAP_RDS <-
    "GSE303006_M23_stage2_EP_mapping_degrees.rds"

MAP_TSV <-
    "GSE303006_M23_stage2_EP_mapping_degrees.tsv.gz"

MANIFEST_RDS <-
    "GSE303006_M23_stage2_implicit_manifest.rds"

HELPER_R <-
    "GSE303006_M23_stage2_implicit_helpers.R"


# ============================================================
# 1. Input checks
# ============================================================

needed <- c(
    EP_FILE,
    STAGE1_FILES
)

missing <- needed[
    !file.exists(needed)
]

if (length(missing) > 0L) {
    stop(
        "Missing input file(s):\n",
        paste(missing, collapse = "\n")
    )
}


# ============================================================
# 2. Load final E-P table
# ============================================================

ep <- as.data.table(
    readRDS(EP_FILE)
)

required_ep_cols <- c(
    "enhancer_id",
    "gene_name",
    "binpair_id"
)

if (!all(required_ep_cols %in% names(ep))) {
    stop(
        "E-P table lacks required column(s): ",
        paste(
            setdiff(
                required_ep_cols,
                names(ep)
            ),
            collapse = ", "
        )
    )
}

if (anyNA(ep$enhancer_id) ||
    anyNA(ep$gene_name) ||
    anyNA(ep$binpair_id)) {
    stop(
        "NA found in enhancer_id, gene_name, or binpair_id."
    )
}

cat(
    "Final E-P pairs :",
    nrow(ep),
    "\n"
)


# ============================================================
# 3. Check E-P uniqueness
#
# enhancer_id x gene_name should be unique in the final set.
# ============================================================

dup_ep <- ep[
    ,
    .N,
    by = .(
        enhancer_id,
        gene_name
    )
][
    N > 1L
]

if (nrow(dup_ep) > 0L) {
    stop(
        "Duplicate enhancer_id x gene_name pairs remain: ",
        nrow(dup_ep)
    )
}


# ============================================================
# 4. Define Stage-2 row ID
#
# Keep original pair_id if it exists, but use a new consecutive
# integer stage2_pair_id for all implicit tensor indexing.
# ============================================================

ep_map <- copy(ep)

ep_map[
    ,
    stage2_pair_id :=
        seq_len(.N)
]

if ("pair_id" %in% names(ep_map)) {

    setcolorder(
        ep_map,
        c(
            "stage2_pair_id",
            "pair_id",
            setdiff(
                names(ep_map),
                c(
                    "stage2_pair_id",
                    "pair_id"
                )
            )
        )
    )

} else {

    setcolorder(
        ep_map,
        c(
            "stage2_pair_id",
            setdiff(
                names(ep_map),
                "stage2_pair_id"
            )
        )
    )
}


# ============================================================
# 5. Compute degrees d_e, d_g, d_b
# ============================================================

deg_e <- ep_map[
    ,
    .(
        d_e = .N
    ),
    by = enhancer_id
]

deg_g <- ep_map[
    ,
    .(
        d_g = .N
    ),
    by = gene_name
]

deg_b <- ep_map[
    ,
    .(
        d_b = .N
    ),
    by = binpair_id
]

ep_map[
    deg_e,
    on = "enhancer_id",
    d_e := i.d_e
]

ep_map[
    deg_g,
    on = "gene_name",
    d_g := i.d_g
]

ep_map[
    deg_b,
    on = "binpair_id",
    d_b := i.d_b
]

stopifnot(
    !anyNA(ep_map$d_e),
    !anyNA(ep_map$d_g),
    !anyNA(ep_map$d_b)
)

stopifnot(
    all(ep_map$d_e >= 1L),
    all(ep_map$d_g >= 1L),
    all(ep_map$d_b >= 1L)
)


# ============================================================
# 6. Degree-correction weights
# ============================================================

ep_map[
    ,
    w_enhancer :=
        1 /
        sqrt(d_e)
]

ep_map[
    ,
    w_gene :=
        1 /
        sqrt(d_g)
]

ep_map[
    ,
    w_binpair :=
        1 /
        sqrt(d_b)
]


# ============================================================
# 7. Load Stage-1 sign-aligned objects
# ============================================================

stage1 <- lapply(
    STAGE1_FILES,
    readRDS
)

names(stage1) <-
    names(STAGE1_FILES)

required_stage1_fields <- c(
    "feature_scores",
    "cell_loadings"
)

for (nm in names(stage1)) {

    obj <- stage1[[nm]]

    if (!all(
        required_stage1_fields %in%
            names(obj)
    )) {
        stop(
            nm,
            " Stage-1 object lacks required field(s): ",
            paste(
                setdiff(
                    required_stage1_fields,
                    names(obj)
                ),
                collapse = ", "
            )
        )
    }

    if (ncol(obj$feature_scores) != K ||
        ncol(obj$cell_loadings) != K) {
        stop(
            nm,
            ": expected K=",
            K,
            "."
        )
    }
}


# ============================================================
# 8. Validate cell order across modalities
# ============================================================

cells <- rownames(
    stage1$ATAC$cell_loadings
)

if (is.null(cells)) {
    stop(
        "ATAC cell_loadings have no row names."
    )
}

for (nm in names(stage1)) {

    current_cells <- rownames(
        stage1[[nm]]$cell_loadings
    )

    if (!identical(
        current_cells,
        cells
    )) {
        stop(
            "Cell ordering differs for modality: ",
            nm
        )
    }
}

n_cell <- length(cells)

cat(
    "Stage-2 cells :",
    n_cell,
    "\n"
)


# ============================================================
# 9. Create integer mapping from each E-P row to Stage-1
#    feature-score row
# ============================================================

B_A <- stage1$ATAC$feature_scores
B_H <- stage1$H3K27ac$feature_scores
B_R <- stage1$RNA$feature_scores
B_D <- stage1$`3D`$feature_scores

if (is.null(rownames(B_A)) ||
    is.null(rownames(B_H)) ||
    is.null(rownames(B_R)) ||
    is.null(rownames(B_D))) {
    stop(
        "One or more feature-score matrices have no row names."
    )
}

# ATAC and H3K27ac both map by enhancer_id.
ep_map[
    ,
    idx_ATAC :=
        match(
            enhancer_id,
            rownames(B_A)
        )
]

ep_map[
    ,
    idx_H3K27ac :=
        match(
            enhancer_id,
            rownames(B_H)
        )
]

# RNA maps by gene_name.
ep_map[
    ,
    idx_RNA :=
        match(
            gene_name,
            rownames(B_R)
        )
]

# 3D maps by binpair_id.
ep_map[
    ,
    idx_3D :=
        match(
            as.character(binpair_id),
            rownames(B_D)
        )
]

stopifnot(
    !anyNA(ep_map$idx_ATAC),
    !anyNA(ep_map$idx_H3K27ac),
    !anyNA(ep_map$idx_RNA),
    !anyNA(ep_map$idx_3D)
)


# ============================================================
# 10. Compact mapping table
#
# Keep the biologically relevant IDs, degrees, weights,
# and integer indices needed by the implicit tensor.
# ============================================================

keep_cols <- c(
    "stage2_pair_id",
    intersect(
        "pair_id",
        names(ep_map)
    ),
    "enhancer_id",
    "gene_name",
    "binpair_id",
    "d_e",
    "d_g",
    "d_b",
    "w_enhancer",
    "w_gene",
    "w_binpair",
    "idx_ATAC",
    "idx_H3K27ac",
    "idx_RNA",
    "idx_3D"
)

mapping <- ep_map[
    ,
    ..keep_cols
]


# ============================================================
# 11. Degree diagnostics
# ============================================================

cat(
    "\nEnhancer degree d_e:\n"
)

print(
    summary(
        deg_e$d_e
    )
)

cat(
    "\nGene degree d_g:\n"
)

print(
    summary(
        deg_g$d_g
    )
)

cat(
    "\n20-kb bin-pair degree d_b:\n"
)

print(
    summary(
        deg_b$d_b
    )
)

cat(
    "\nMaximum degrees:\n"
)

cat(
    "max d_e =",
    max(deg_e$d_e),
    "\n"
)

cat(
    "max d_g =",
    max(deg_g$d_g),
    "\n"
)

cat(
    "max d_b =",
    max(deg_b$d_b),
    "\n"
)


# ============================================================
# 12. Verify the degree-correction identity
#
# For every underlying enhancer e:
#   sum_{p containing e} w_enhancer[p]^2 = 1
#
# Similarly for genes and bin-pairs.
# ============================================================

check_e <- mapping[
    ,
    .(
        sum_w2 =
            sum(
                w_enhancer^2
            )
    ),
    by = enhancer_id
]

check_g <- mapping[
    ,
    .(
        sum_w2 =
            sum(
                w_gene^2
            )
    ),
    by = gene_name
]

check_b <- mapping[
    ,
    .(
        sum_w2 =
            sum(
                w_binpair^2
            )
    ),
    by = binpair_id
]

tol <- 1e-10

if (max(abs(check_e$sum_w2 - 1)) > tol) {
    stop(
        "Enhancer degree-correction identity failed."
    )
}

if (max(abs(check_g$sum_w2 - 1)) > tol) {
    stop(
        "Gene degree-correction identity failed."
    )
}

if (max(abs(check_b$sum_w2 - 1)) > tol) {
    stop(
        "Bin-pair degree-correction identity failed."
    )
}

cat(
    "\nDegree-correction identities verified.\n"
)


# ============================================================
# 13. Save mapping
# ============================================================

saveRDS(
    mapping,
    MAP_RDS
)

fwrite(
    mapping,
    MAP_TSV,
    sep = "\t"
)


# ============================================================
# 14. Save compact implicit-tensor manifest
#
# Do not duplicate the Stage-1 matrices into this object.
# Save only file references + dimensions + conventions.
# ============================================================

manifest <- list(

    tensor_name =
        "GSE303006_M23_hierarchical_EP_stage2",

    dimensions = c(
        EP_pair = nrow(mapping),
        cell = n_cell,
        component = K,
        modality = 4L
    ),

    dimension_order = c(
        "EP_pair",
        "cell",
        "component",
        "modality"
    ),

    modalities = c(
        "ATAC",
        "H3K27ac",
        "RNA",
        "3D"
    ),

    K = K,

    cells = cells,

    mapping_file = MAP_RDS,

    stage1_files = STAGE1_FILES,

    degree_correction = list(
        ATAC = "1/sqrt(d_e)",
        H3K27ac = "1/sqrt(d_e)",
        RNA = "1/sqrt(d_g)",
        `3D` = "1/sqrt(d_b)"
    ),

    element_definition = c(
        ATAC =
            "w_enhancer[p] * B_ATAC[e(p),k] * V_ATAC[c,k]",
        H3K27ac =
            "w_enhancer[p] * B_H3K27ac[e(p),k] * V_H3K27ac[c,k]",
        RNA =
            "w_gene[p] * B_RNA[g(p),k] * V_RNA[c,k]",
        `3D` =
            "w_binpair[p] * B_3D[b(p),k] * V_3D[c,k]"
    ),

    materialized = FALSE
)

saveRDS(
    manifest,
    MANIFEST_RDS
)


# ============================================================
# 15. Write helper functions as a standalone R file
# ============================================================

helper_code <- paste0(
'
# ============================================================
# Implicit Stage-2 tensor helper functions
# Auto-generated by 07_prepare_stage2_implicit.R
# ============================================================

library(data.table)

load_stage2_implicit <- function(
    manifest_file =
        "GSE303006_M23_stage2_implicit_manifest.rds"
) {

    M <- readRDS(
        manifest_file
    )

    map <- as.data.table(
        readRDS(
            M$mapping_file
        )
    )

    S <- lapply(
        M$stage1_files,
        readRDS
    )

    names(S) <-
        names(
            M$stage1_files
        )

    # Validate cells.
    for (nm in names(S)) {

        if (!identical(
            rownames(
                S[[nm]]$cell_loadings
            ),
            M$cells
        )) {
            stop(
                "Cell order mismatch in ",
                nm
            )
        }
    }

    list(
        manifest = M,
        map = map,
        stage1 = S
    )
}


stage2_dims <- function(S) {

    S$manifest$dimensions
}


stage2_block <- function(
    S,
    pair_idx,
    cell_idx,
    k_idx = seq_len(
        S$manifest$K
    ),
    modalities =
        S$manifest$modalities
) {

    pair_idx <- as.integer(
        pair_idx
    )

    cell_idx <- as.integer(
        cell_idx
    )

    k_idx <- as.integer(
        k_idx
    )

    if (any(
        pair_idx < 1L |
        pair_idx > nrow(S$map)
    )) {
        stop(
            "pair_idx out of range."
        )
    }

    if (any(
        cell_idx < 1L |
        cell_idx >
            length(
                S$manifest$cells
            )
    )) {
        stop(
            "cell_idx out of range."
        )
    }

    if (any(
        k_idx < 1L |
        k_idx >
            S$manifest$K
    )) {
        stop(
            "k_idx out of range."
        )
    }

    bad_mod <- setdiff(
        modalities,
        S$manifest$modalities
    )

    if (length(bad_mod) > 0L) {
        stop(
            "Unknown modality: ",
            paste(
                bad_mod,
                collapse = ", "
            )
        )
    }

    mp <- S$map[
        pair_idx
    ]

    out <- array(
        0,
        dim = c(
            length(pair_idx),
            length(cell_idx),
            length(k_idx),
            length(modalities)
        ),
        dimnames = list(
            EP_pair =
                as.character(
                    mp$stage2_pair_id
                ),
            cell =
                S$manifest$cells[
                    cell_idx
                ],
            component =
                sprintf(
                    "K%02d",
                    k_idx
                ),
            modality =
                modalities
        )
    )

    for (mm in seq_along(
        modalities
    )) {

        mod <- modalities[mm]

        if (mod == "ATAC") {

            B <- S$stage1$ATAC$
                feature_scores[
                    mp$idx_ATAC,
                    k_idx,
                    drop = FALSE
                ]

            V <- S$stage1$ATAC$
                cell_loadings[
                    cell_idx,
                    k_idx,
                    drop = FALSE
                ]

            W <- mp$w_enhancer

        } else if (
            mod == "H3K27ac"
        ) {

            B <- S$stage1$H3K27ac$
                feature_scores[
                    mp$idx_H3K27ac,
                    k_idx,
                    drop = FALSE
                ]

            V <- S$stage1$H3K27ac$
                cell_loadings[
                    cell_idx,
                    k_idx,
                    drop = FALSE
                ]

            W <- mp$w_enhancer

        } else if (
            mod == "RNA"
        ) {

            B <- S$stage1$RNA$
                feature_scores[
                    mp$idx_RNA,
                    k_idx,
                    drop = FALSE
                ]

            V <- S$stage1$RNA$
                cell_loadings[
                    cell_idx,
                    k_idx,
                    drop = FALSE
                ]

            W <- mp$w_gene

        } else if (
            mod == "3D"
        ) {

            B <- S$stage1$`3D`$
                feature_scores[
                    mp$idx_3D,
                    k_idx,
                    drop = FALSE
                ]

            V <- S$stage1$`3D`$
                cell_loadings[
                    cell_idx,
                    k_idx,
                    drop = FALSE
                ]

            W <- mp$w_binpair

        } else {

            stop(
                "Internal modality error."
            )
        }

        # For each k, tensor slice is
        #   (W * B[,k]) outer V[,k]
        for (jj in seq_along(
            k_idx
        )) {

            out[
                ,
                ,
                jj,
                mm
            ] <- tcrossprod(
                W *
                    B[
                        ,
                        jj
                    ],
                V[
                    ,
                    jj
                ]
            )
        }
    }

    out
}


# Iterate through pair blocks without materializing all X.
stage2_pair_blocks <- function(
    S,
    block_size = 1000L
) {

    n <- nrow(
        S$map
    )

    split(
        seq_len(n),
        ceiling(
            seq_len(n) /
            as.integer(
                block_size
            )
        )
    )
}


# Example:
#
#   S <- load_stage2_implicit()
#   stage2_dims(S)
#
#   Xsmall <- stage2_block(
#       S,
#       pair_idx = 1:10,
#       cell_idx = 1:20,
#       k_idx = 1:3,
#       modalities = c(
#           "ATAC",
#           "RNA"
#       )
#   )
#
#   dim(Xsmall)
#
# This constructs only:
#   10 x 20 x 3 x 2
# not the full tensor.
'
)

writeLines(
    helper_code,
    HELPER_R
)


# ============================================================
# 16. Small read-back / block-generation test
# ============================================================

source(
    HELPER_R
)

S_test <- load_stage2_implicit(
    MANIFEST_RDS
)

cat(
    "\nImplicit tensor dimensions:\n"
)

print(
    stage2_dims(
        S_test
    )
)

X_test <- stage2_block(
    S_test,
    pair_idx = 1:3,
    cell_idx = 1:5,
    k_idx = 1:2,
    modalities = c(
        "ATAC",
        "H3K27ac",
        "RNA",
        "3D"
    )
)

cat(
    "\nSmall generated block dimension:\n"
)

print(
    dim(
        X_test
    )
)

cat(
    "\nSmall generated block summary:\n"
)

print(
    summary(
        as.numeric(
            X_test
        )
    )
)


# ============================================================
# 17. Final summary
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Stage-2 implicit representation prepared\n"
)

cat(
    "============================================================\n"
)

cat(
    "E-P pairs    :",
    nrow(mapping),
    "\n"
)

cat(
    "Cells        :",
    n_cell,
    "\n"
)

cat(
    "Components   :",
    K,
    "\n"
)

cat(
    "Modalities   : 4\n"
)

cat(
    "\nConceptual tensor:\n"
)

cat(
    nrow(mapping),
    " x ",
    n_cell,
    " x ",
    K,
    " x 4\n",
    sep = ""
)

cat(
    "\nThe full tensor was NOT materialized.\n"
)

cat(
    "Use ",
    HELPER_R,
    " and stage2_block() to generate only required blocks.\n",
    sep = ""
)

cat(
    "\nSaved mapping:\n",
    MAP_RDS,
    "\n",
    MAP_TSV,
    "\n",
    sep = ""
)

cat(
    "\nSaved manifest:\n",
    MANIFEST_RDS,
    "\n",
    sep = ""
)
