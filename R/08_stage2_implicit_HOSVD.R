# ============================================================
# 08_stage2_implicit_tucker.R
# GSE303006 / GENCODE M23
#
# Stage 2 Tucker/HOSVD on the implicit four-way tensor:
#
#   X[p, c, k, m]
#   730,969 E-P x 4,258 cells x 20 Stage-1 components x 4 modalities
#
# WITHOUT materializing X.
#
# Input element definition:
#
#   X[p,c,k,ATAC] =
#     w_e[p] * B_A[e(p),k] * V_A[c,k]
#
#   X[p,c,k,H3K27ac] =
#     w_e[p] * B_H[e(p),k] * V_H[c,k]
#
#   X[p,c,k,RNA] =
#     w_g[p] * B_R[g(p),k] * V_R[c,k]
#
#   X[p,c,k,3D] =
#     w_b[p] * B_D[b(p),k] * V_D[c,k]
#
# Key fact:
#   The tensor lies in at most an 80-dimensional E-P subspace
#   and an 80-dimensional cell subspace because there are only
#   20 rank-1 Stage-1 components x 4 modalities.
#
# Therefore:
#   - no 730,969 x 730,969 Gram matrix
#   - no full 730,969 x 4,258 x 20 x 4 tensor
#
# Default Tucker ranks:
#   E-P       : 20
#   cell      : 20
#   component : 20 (kept full; no mixing)
#   modality  : 4  (full HOSVD rotation)
#
# These are settings and can be changed later.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
}

library(data.table)


# ============================================================
# 1. Settings
# ============================================================

R_EP <- 20L
R_CELL <- 20L
R_K <- 20L
R_MODAL <- 4L

PAIR_BLOCK <- 10000L

MANIFEST_FILE <-
    "GSE303006_M23_stage2_implicit_manifest.rds"

MAPPING_FILE <-
    "GSE303006_M23_stage2_EP_mapping_degrees.rds"

OUT_EP_FACTOR <-
    "GSE303006_M23_stage2_Tucker_EP_factor.rds"

OUT_SMALL_FACTORS <-
    "GSE303006_M23_stage2_Tucker_small_factors.rds"

OUT_CORE <-
    "GSE303006_M23_stage2_Tucker_core.rds"

OUT_DIAG <-
    "GSE303006_M23_stage2_Tucker_diagnostics.rds"

OUT_MODAL_TSV <-
    "GSE303006_M23_stage2_Tucker_modality_factor.tsv"

OUT_CELL_TSV <-
    "GSE303006_M23_stage2_Tucker_cell_factor.tsv"


# ============================================================
# 2. Load manifest / mapping / Stage-1 objects
# ============================================================

if (!file.exists(MANIFEST_FILE)) {
    stop("Missing manifest: ", MANIFEST_FILE)
}

if (!file.exists(MAPPING_FILE)) {
    stop("Missing mapping: ", MAPPING_FILE)
}

M <- readRDS(
    MANIFEST_FILE
)

map <- as.data.table(
    readRDS(
        MAPPING_FILE
    )
)

S <- lapply(
    M$stage1_files,
    readRDS
)

names(S) <- names(
    M$stage1_files
)

required_modalities <- c(
    "ATAC",
    "H3K27ac",
    "RNA",
    "3D"
)

if (!identical(
    names(S),
    required_modalities
)) {
    # reorder if all names are present
    if (!all(
        required_modalities %in%
            names(S)
    )) {
        stop(
            "Stage-1 modalities do not match expected set."
        )
    }

    S <- S[
        required_modalities
    ]
}

K <- as.integer(
    M$K
)

cells <- as.character(
    M$cells
)

P <- nrow(
    map
)

C <- length(
    cells
)

MOD <- length(
    required_modalities
)

stopifnot(
    K == 20L,
    MOD == 4L,
    R_EP >= 1L,
    R_EP <= 4L * K,
    R_CELL >= 1L,
    R_CELL <= 4L * K,
    R_K == K,
    R_MODAL >= 1L,
    R_MODAL <= MOD
)

cat(
    "Implicit tensor dimensions: ",
    P, " x ", C, " x ", K, " x ", MOD,
    "\n",
    sep = ""
)

cat(
    "Requested Tucker ranks: ",
    R_EP, " x ", R_CELL, " x ", R_K, " x ", R_MODAL,
    "\n",
    sep = ""
)


# ============================================================
# 3. Validate Stage-1 objects
# ============================================================

for (nm in required_modalities) {

    obj <- S[[nm]]

    if (!all(
        c(
            "feature_scores",
            "cell_loadings"
        ) %in%
            names(obj)
    )) {
        stop(
            nm,
            " lacks feature_scores/cell_loadings."
        )
    }

    if (ncol(
        obj$feature_scores
    ) != K) {
        stop(
            nm,
            " feature_scores K mismatch."
        )
    }

    if (ncol(
        obj$cell_loadings
    ) != K) {
        stop(
            nm,
            " cell_loadings K mismatch."
        )
    }

    if (!identical(
        rownames(
            obj$cell_loadings
        ),
        cells
    )) {
        stop(
            nm,
            " cell ordering mismatch."
        )
    }

    ortho_err <- max(
        abs(
            crossprod(
                obj$cell_loadings
            ) -
            diag(K)
        )
    )

    cat(
        nm,
        " cell orthogonality error = ",
        format(
            ortho_err,
            scientific = TRUE
        ),
        "\n",
        sep = ""
    )
}


# ============================================================
# 4. Aliases
# ============================================================

BA <- S$ATAC$feature_scores
BH <- S$H3K27ac$feature_scores
BR <- S$RNA$feature_scores
BD <- S$`3D`$feature_scores

VA <- S$ATAC$cell_loadings
VH <- S$H3K27ac$cell_loadings
VR <- S$RNA$cell_loadings
VD <- S$`3D`$cell_loadings

stopifnot(
    !anyNA(map$idx_ATAC),
    !anyNA(map$idx_H3K27ac),
    !anyNA(map$idx_RNA),
    !anyNA(map$idx_3D)
)


# ============================================================
# 5. Helper: build only the pair x 80 Stage-1 amplitude block
#
# Columns:
#   1:20   ATAC
#   21:40  H3K27ac
#   41:60  RNA
#   61:80  3D
#
# A[p,(m,k)] is the pair-side amplitude after degree correction.
# ============================================================

build_A_block <- function(
    ii
) {

    mp <- map[
        ii
    ]

    A1 <- BA[
        mp$idx_ATAC,
        ,
        drop = FALSE
    ]

    A1 <- A1 *
        mp$w_enhancer

    A2 <- BH[
        mp$idx_H3K27ac,
        ,
        drop = FALSE
    ]

    A2 <- A2 *
        mp$w_enhancer

    A3 <- BR[
        mp$idx_RNA,
        ,
        drop = FALSE
    ]

    A3 <- A3 *
        mp$w_gene

    A4 <- BD[
        mp$idx_3D,
        ,
        drop = FALSE
    ]

    A4 <- A4 *
        mp$w_binpair

    cbind(
        A1,
        A2,
        A3,
        A4
    )
}


# ============================================================
# 6. Pair-side 80 x 80 Gram matrix
#
# H_A = A^T A.
#
# Non-zero eigenvalues of H_A are exactly the non-zero
# eigenvalues of the enormous E-P mode Gram matrix A A^T.
# ============================================================

L <- K * MOD

H_A <- matrix(
    0,
    nrow = L,
    ncol = L
)

pair_blocks <- split(
    seq_len(P),
    ceiling(
        seq_len(P) /
        PAIR_BLOCK
    )
)

cat(
    "\nComputing exact 80 x 80 pair-side Gram matrix in ",
    length(pair_blocks),
    " blocks...\n",
    sep = ""
)

for (ib in seq_along(
    pair_blocks
)) {

    ii <- pair_blocks[[ib]]

    Ab <- build_A_block(
        ii
    )

    H_A <- H_A +
        crossprod(
            Ab
        )

    if (
        ib %% 10L == 0L ||
        ib == length(pair_blocks)
    ) {
        cat(
            "pair Gram block ",
            ib,
            "/",
            length(pair_blocks),
            "\n",
            sep = ""
        )
    }

    rm(
        Ab
    )
}

H_A <- (
    H_A +
    t(H_A)
) / 2


# ============================================================
# 7. Tensor total norm and modality norms
#
# With Stage-1 retained norms normalized to 1 and the
# degree-correction identities, each modality should contribute
# Frobenius norm^2 ~= 1, hence total tensor norm^2 ~= 4.
# ============================================================

offset <- list(
    ATAC = 1:K,
    H3K27ac = (K + 1):(2 * K),
    RNA = (2 * K + 1):(3 * K),
    `3D` = (3 * K + 1):(4 * K)
)

modal_norm2 <- sapply(
    offset,
    function(jj) {
        sum(
            diag(
                H_A[
                    jj,
                    jj,
                    drop = FALSE
                ]
            )
        )
    }
)

tensor_fro2 <- sum(
    modal_norm2
)

cat(
    "\nStage-2 modality Frobenius norm^2:\n"
)

print(
    modal_norm2
)

cat(
    "Total tensor Frobenius norm^2 = ",
    tensor_fro2,
    "\n",
    sep = ""
)


# ============================================================
# 8. E-P mode factor from the 80 x 80 Gram matrix
# ============================================================

eeP <- eigen(
    H_A,
    symmetric = TRUE
)

lambdaP_all <- pmax(
    as.numeric(
        eeP$values
    ),
    0
)

WP_all <- eeP$vectors

positiveP <- which(
    lambdaP_all >
        max(lambdaP_all) *
        1e-12
)

effective_rank_P <- length(
    positiveP
)

if (R_EP > effective_rank_P) {
    stop(
        "R_EP exceeds numerical E-P mode rank."
    )
}

lambdaP <- lambdaP_all[
    seq_len(R_EP)
]

sigmaP <- sqrt(
    lambdaP
)

WP <- WP_all[
    ,
    seq_len(R_EP),
    drop = FALSE
]

cat(
    "\nE-P mode numerical rank <= ",
    effective_rank_P,
    "\n",
    sep = ""
)

cat(
    "Top E-P mode singular values:\n"
)

print(
    sqrt(
        lambdaP_all[
            seq_len(
                min(
                    80L,
                    length(lambdaP_all)
                )
            )
        ]
    )
)

ep_mode_fraction <- sum(
    lambdaP
) /
    tensor_fro2

cat(
    "E-P mode top-",
    R_EP,
    " energy fraction = ",
    ep_mode_fraction,
    "\n",
    sep = ""
)


# ============================================================
# 9. Build and save E-P factor U_EP blockwise
#
# If A = U Sigma W^T, then
#   U = A W Sigma^{-1}
# ============================================================

U_EP <- matrix(
    0,
    nrow = P,
    ncol = R_EP
)

cat(
    "\nBuilding E-P factor ",
    P, " x ", R_EP,
    " blockwise...\n",
    sep = ""
)

for (ib in seq_along(
    pair_blocks
)) {

    ii <- pair_blocks[[ib]]

    Ab <- build_A_block(
        ii
    )

    Ub <- Ab %*%
        WP

    Ub <- sweep(
        Ub,
        2,
        sigmaP,
        "/"
    )

    U_EP[
        ii,
    ] <- Ub

    if (
        ib %% 10L == 0L ||
        ib == length(pair_blocks)
    ) {
        cat(
            "E-P factor block ",
            ib,
            "/",
            length(pair_blocks),
            "\n",
            sep = ""
        )
    }

    rm(
        Ab,
        Ub
    )
}

rownames(
    U_EP
) <- as.character(
    map$stage2_pair_id
)

colnames(
    U_EP
) <- sprintf(
    "EP%02d",
    seq_len(R_EP)
)

ep_orth_err <- max(
    abs(
        crossprod(
            U_EP
        ) -
        diag(
            R_EP
        )
    )
)

cat(
    "E-P factor orthogonality error = ",
    format(
        ep_orth_err,
        scientific = TRUE
    ),
    "\n",
    sep = ""
)

saveRDS(
    U_EP,
    OUT_EP_FACTOR
)


# ============================================================
# 10. Cell-mode factor
#
# Cell mode Gram:
#
#   G_C =
#     sum_{m,k} ||a_{m,k}||^2 v_{m,k} v_{m,k}^T
#
# Construct Ccat (4258 x 80), whose columns are
#   ||a_{m,k}|| * v_{m,k}
#
# Then U_CELL is the left singular basis of Ccat.
# ============================================================

normA <- sqrt(
    pmax(
        diag(
            H_A
        ),
        0
    )
)

Ccat <- cbind(
    sweep(
        VA,
        2,
        normA[
            offset$ATAC
        ],
        "*"
    ),
    sweep(
        VH,
        2,
        normA[
            offset$H3K27ac
        ],
        "*"
    ),
    sweep(
        VR,
        2,
        normA[
            offset$RNA
        ],
        "*"
    ),
    sweep(
        VD,
        2,
        normA[
            offset$`3D`
        ],
        "*"
    )
)

svC <- svd(
    Ccat,
    nu = R_CELL,
    nv = 0
)

U_CELL <- svC$u[
    ,
    seq_len(R_CELL),
    drop = FALSE
]

sigmaC_all <- svC$d

sigmaC <- sigmaC_all[
    seq_len(R_CELL)
]

rownames(
    U_CELL
) <- cells

colnames(
    U_CELL
) <- sprintf(
    "CELL%02d",
    seq_len(R_CELL)
)

cell_mode_fraction <- sum(
    sigmaC^2
) /
    tensor_fro2

cat(
    "\nTop cell-mode singular values:\n"
)

print(
    sigmaC_all
)

cat(
    "Cell mode top-",
    R_CELL,
    " energy fraction = ",
    cell_mode_fraction,
    "\n",
    sep = ""
)


# ============================================================
# 11. Component (K) mode
#
# Because each Stage-1 V_m has orthonormal columns,
# the K-mode Gram is diagonal:
#
#   G_K[k,k'] = 0 for k != k'
#
# We deliberately KEEP the original K=1..20 ordering and use
# identity U_K. This avoids unnecessary permutation/mixing.
# ============================================================

energy_K <- numeric(
    K
)

for (k in seq_len(K)) {

    energy_K[k] <-
        H_A[
            offset$ATAC[k],
            offset$ATAC[k]
        ] +
        H_A[
            offset$H3K27ac[k],
            offset$H3K27ac[k]
        ] +
        H_A[
            offset$RNA[k],
            offset$RNA[k]
        ] +
        H_A[
            offset$`3D`[k],
            offset$`3D`[k]
        ]
}

U_K <- diag(
    K
)

rownames(
    U_K
) <- sprintf(
    "K%02d",
    seq_len(K)
)

colnames(
    U_K
) <- rownames(
    U_K
)

cat(
    "\nComponent-mode energies (original K order preserved):\n"
)

print(
    energy_K
)


# ============================================================
# 12. Modality-mode Gram matrix and factor
#
# G_M[m,n] =
#   sum_k
#     <a_mk, a_nk> *
#     <v_mk, v_nk>
#
# ============================================================

Vlist <- list(
    ATAC = VA,
    H3K27ac = VH,
    RNA = VR,
    `3D` = VD
)

G_M <- matrix(
    0,
    nrow = MOD,
    ncol = MOD,
    dimnames = list(
        required_modalities,
        required_modalities
    )
)

for (im in seq_len(
    MOD
)) {

    m1 <- required_modalities[
        im
    ]

    j1 <- offset[[m1]]

    for (jm in seq_len(
        MOD
    )) {

        m2 <- required_modalities[
            jm
        ]

        j2 <- offset[[m2]]

        pair_cross <-
            diag(
                H_A[
                    j1,
                    j2,
                    drop = FALSE
                ]
            )

        cell_cross <-
            diag(
                crossprod(
                    Vlist[[m1]],
                    Vlist[[m2]]
                )
            )

        G_M[
            im,
            jm
        ] <- sum(
            pair_cross *
                cell_cross
        )
    }
}

G_M <- (
    G_M +
    t(
        G_M
    )
) / 2

eeM <- eigen(
    G_M,
    symmetric = TRUE
)

lambdaM_all <- pmax(
    as.numeric(
        eeM$values
    ),
    0
)

U_MODAL <- eeM$vectors[
    ,
    seq_len(R_MODAL),
    drop = FALSE
]

rownames(
    U_MODAL
) <- required_modalities

colnames(
    U_MODAL
) <- sprintf(
    "MOD%02d",
    seq_len(R_MODAL)
)

# deterministic sign
for (r in seq_len(
    R_MODAL
)) {

    jj <- which.max(
        abs(
            U_MODAL[
                ,
                r
            ]
        )
    )

    if (
        U_MODAL[
            jj,
            r
        ] < 0
    ) {
        U_MODAL[
            ,
            r
        ] <- -
            U_MODAL[
                ,
                r
            ]
    }
}

modal_mode_fraction <- sum(
    lambdaM_all[
        seq_len(R_MODAL)
    ]
) /
    tensor_fro2

cat(
    "\nModality-mode Gram matrix:\n"
)

print(
    G_M
)

cat(
    "\nModality-mode eigenvalues:\n"
)

print(
    lambdaM_all
)

cat(
    "\nModality factor:\n"
)

print(
    U_MODAL
)

cat(
    "Modality mode top-",
    R_MODAL,
    " energy fraction = ",
    modal_mode_fraction,
    "\n",
    sep = ""
)


# ============================================================
# 13. Core tensor WITHOUT materializing X
#
# Pair projection:
#   U_EP^T A = Sigma_P W_P^T
#
# This is R_EP x 80.
# ============================================================

PAIR_PROJ <- diag(
    sigmaP,
    nrow = R_EP,
    ncol = R_EP
) %*%
    t(
        WP
    )

# Cell projections U_CELL^T V_m
CELL_PROJ <- lapply(
    Vlist,
    function(Vm) {
        crossprod(
            U_CELL,
            Vm
        )
    }
)

# Pre-modality core:
#   R_EP x R_CELL x K x 4
#
# K mode is identity / original order.
core_pre <- array(
    0,
    dim = c(
        R_EP,
        R_CELL,
        K,
        MOD
    ),
    dimnames = list(
        EP = colnames(U_EP),
        CELL = colnames(U_CELL),
        K = sprintf(
            "K%02d",
            seq_len(K)
        ),
        MODALITY = required_modalities
    )
)

for (im in seq_len(
    MOD
)) {

    nm <- required_modalities[
        im
    ]

    jj <- offset[[nm]]

    Cp <- CELL_PROJ[[nm]]

    for (k in seq_len(
        K
    )) {

        pair_vec <-
            PAIR_PROJ[
                ,
                jj[k]
            ]

        cell_vec <-
            Cp[
                ,
                k
            ]

        core_pre[
            ,
            ,
            k,
            im
        ] <- tcrossprod(
            pair_vec,
            cell_vec
        )
    }
}

# Apply modality factor U_MODAL.
CORE <- array(
    0,
    dim = c(
        R_EP,
        R_CELL,
        R_K,
        R_MODAL
    ),
    dimnames = list(
        EP = colnames(U_EP),
        CELL = colnames(U_CELL),
        K = sprintf(
            "K%02d",
            seq_len(R_K)
        ),
        MODAL = colnames(U_MODAL)
    )
)

for (rmod in seq_len(
    R_MODAL
)) {

    for (im in seq_len(
        MOD
    )) {

        CORE[
            ,
            ,
            ,
            rmod
        ] <- CORE[
            ,
            ,
            ,
            rmod
        ] +
            core_pre[
                ,
                ,
                seq_len(R_K),
                im
            ] *
            U_MODAL[
                im,
                rmod
            ]
    }
}

core_fro2 <- sum(
    CORE^2
)

multilinear_fraction <-
    core_fro2 /
    tensor_fro2

cat(
    "\nCore Frobenius norm^2 = ",
    core_fro2,
    "\n",
    sep = ""
)

cat(
    "Multilinear Tucker energy fraction = ",
    multilinear_fraction,
    "\n",
    sep = ""
)


# ============================================================
# 14. Save small factors / core / diagnostics
# ============================================================

small_factors <- list(

    ranks = c(
        EP = R_EP,
        CELL = R_CELL,
        K = R_K,
        MODALITY = R_MODAL
    ),

    U_CELL = U_CELL,

    U_K = U_K,

    U_MODAL = U_MODAL,

    cells = cells,

    modalities = required_modalities,

    component_labels =
        sprintf(
            "K%02d",
            seq_len(K)
        )
)

saveRDS(
    small_factors,
    OUT_SMALL_FACTORS
)

saveRDS(
    CORE,
    OUT_CORE
)

diagnostics <- list(

    tensor_dimensions = c(
        EP = P,
        CELL = C,
        K = K,
        MODALITY = MOD
    ),

    tucker_ranks = c(
        EP = R_EP,
        CELL = R_CELL,
        K = R_K,
        MODALITY = R_MODAL
    ),

    tensor_frobenius2 =
        tensor_fro2,

    modality_frobenius2 =
        modal_norm2,

    pair_mode_singular_values =
        sqrt(
            lambdaP_all
        ),

    cell_mode_singular_values =
        sigmaC_all,

    component_mode_energy =
        energy_K,

    modality_mode_gram =
        G_M,

    modality_mode_eigenvalues =
        lambdaM_all,

    EP_mode_energy_fraction =
        ep_mode_fraction,

    cell_mode_energy_fraction =
        cell_mode_fraction,

    modality_mode_energy_fraction =
        modal_mode_fraction,

    core_frobenius2 =
        core_fro2,

    multilinear_energy_fraction =
        multilinear_fraction,

    EP_factor_orthogonality_error =
        ep_orth_err
)

saveRDS(
    diagnostics,
    OUT_DIAG
)

fwrite(
    as.data.table(
        U_MODAL,
        keep.rownames =
            "modality"
    ),
    OUT_MODAL_TSV,
    sep = "\t"
)

fwrite(
    as.data.table(
        U_CELL,
        keep.rownames =
            "cell"
    ),
    OUT_CELL_TSV,
    sep = "\t"
)


# ============================================================
# 15. Final summary
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Stage-2 implicit Tucker/HOSVD completed\n"
)

cat(
    "============================================================\n"
)

cat(
    "Original conceptual tensor: ",
    P, " x ", C, " x ", K, " x ", MOD,
    "\n",
    sep = ""
)

cat(
    "Tucker core: ",
    paste(
        dim(CORE),
        collapse = " x "
    ),
    "\n",
    sep = ""
)

cat(
    "E-P factor: ",
    paste(
        dim(U_EP),
        collapse = " x "
    ),
    "\n",
    sep = ""
)

cat(
    "Cell factor: ",
    paste(
        dim(U_CELL),
        collapse = " x "
    ),
    "\n",
    sep = ""
)

cat(
    "Component factor: ",
    paste(
        dim(U_K),
        collapse = " x "
    ),
    " (identity; original K order retained)\n",
    sep = ""
)

cat(
    "Modality factor: ",
    paste(
        dim(U_MODAL),
        collapse = " x "
    ),
    "\n",
    sep = ""
)

cat(
    "\nEnergy captured by selected multilinear subspace: ",
    sprintf(
        "%.4f",
        multilinear_fraction
    ),
    "\n",
    sep = ""
)

cat(
    "\nSaved files:\n",
    OUT_EP_FACTOR, "\n",
    OUT_SMALL_FACTORS, "\n",
    OUT_CORE, "\n",
    OUT_DIAG, "\n",
    OUT_MODAL_TSV, "\n",
    OUT_CELL_TSV, "\n",
    sep = ""
)
