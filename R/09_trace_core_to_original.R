# ============================================================
# 10_trace_core_to_original.R
# GSE303006 / GENCODE M23
#
# Trace the NON-K-MIXED Stage-2 Tucker core back to:
#   1) original Stage-1 (modality, K) directions
#   2) original modalities
#   3) original E-P pairs
#   4) original cells
#
# Uses ONLY the K-unmixed output from 08_stage2_implicit_tucker.R.
#
# Core:
#   CORE[a,b,k,r]
#
# where
#   a = E-P HOSVD component
#   b = cell HOSVD component
#   k = original Stage-1 K1...K20 (NOT mixed)
#   r = modality HOSVD component
#
# Key reverse projections:
#
#   PairProj[a,m,k] =
#       < U_EP[,a], A_{m,k} >
#
#   CellProj[b,m,k] =
#       < U_CELL[,b], V_{m,k} >
#
#   CorePre[a,b,k,m] =
#       PairProj[a,m,k] * CellProj[b,m,k]
#
#   CORE[a,b,k,r] =
#       sum_m CorePre[a,b,k,m] * U_MODAL[m,r]
#
# Thus, for every large CORE[a,b,k,r], we can quantify
# exactly how ATAC / H3K27ac / RNA / 3D contributed to it.
#
# We also rank the original E-P pairs and cells associated with
# the selected HOSVD components.
# ============================================================


# ============================================================
# 0. Packages / settings
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
}

library(data.table)

TOP_CORE <- 100L
TOP_EP_EACH_SIGN <- 50L
TOP_CELL_EACH_SIGN <- 50L
TOP_ORIGIN_80 <- 10L
PAIR_BLOCK <- 10000L

CORE_FILE <-
    "GSE303006_M23_stage2_Tucker_core.rds"

SMALL_FACTORS_FILE <-
    "GSE303006_M23_stage2_Tucker_small_factors.rds"

EP_FACTOR_FILE <-
    "GSE303006_M23_stage2_Tucker_EP_factor.rds"

MAPPING_FILE <-
    "GSE303006_M23_stage2_EP_mapping_degrees.rds"

MANIFEST_FILE <-
    "GSE303006_M23_stage2_implicit_manifest.rds"

OUT_TOP_CORE <-
    "GSE303006_M23_stage2_top100_core.tsv"

OUT_MODAL_DECOMP <-
    "GSE303006_M23_stage2_top100_core_original_modality_decomposition.tsv"

OUT_EP_ORIGIN <-
    "GSE303006_M23_stage2_selected_EPcomponents_origin_80.tsv"

OUT_CELL_ORIGIN <-
    "GSE303006_M23_stage2_selected_CELLcomponents_origin_80.tsv"

OUT_TOP_EP <-
    "GSE303006_M23_stage2_top100_core_top_EP.tsv.gz"

OUT_TOP_CELL <-
    "GSE303006_M23_stage2_top100_core_top_cells.tsv.gz"

OUT_CORE_PRE <-
    "GSE303006_M23_stage2_core_before_modality_rotation.rds"

OUT_TRACE_RDS <-
    "GSE303006_M23_stage2_core_trace_results.rds"


# ============================================================
# 1. Load
# ============================================================

need <- c(
    CORE_FILE,
    SMALL_FACTORS_FILE,
    EP_FACTOR_FILE,
    MAPPING_FILE,
    MANIFEST_FILE
)

missing <- need[!file.exists(need)]

if (length(missing) > 0L) {
    stop(
        "Missing input file(s):\n",
        paste(missing, collapse = "\n")
    )
}

CORE <- readRDS(CORE_FILE)
F <- readRDS(SMALL_FACTORS_FILE)
U_EP <- readRDS(EP_FACTOR_FILE)
map <- as.data.table(readRDS(MAPPING_FILE))
M <- readRDS(MANIFEST_FILE)

S <- lapply(
    M$stage1_files,
    readRDS
)

names(S) <- names(M$stage1_files)

mods <- c(
    "ATAC",
    "H3K27ac",
    "RNA",
    "3D"
)

if (!all(mods %in% names(S))) {
    stop("Stage-1 modality names do not match.")
}

S <- S[mods]

U_CELL <- F$U_CELL
U_MODAL <- F$U_MODAL

if (is.null(U_CELL) || is.null(U_MODAL)) {
    stop("U_CELL or U_MODAL missing from small_factors.")
}

if (length(dim(CORE)) != 4L) {
    stop("CORE must be 4-way.")
}

R_EP <- dim(CORE)[1]
R_CELL <- dim(CORE)[2]
K <- dim(CORE)[3]
R_MODAL <- dim(CORE)[4]

P <- nrow(map)
C <- nrow(U_CELL)

cat(
    "CORE dimensions :",
    paste(dim(CORE), collapse = " x "),
    "\n"
)

cat(
    "U_EP dimensions :",
    paste(dim(U_EP), collapse = " x "),
    "\n"
)

cat(
    "U_CELL dimensions :",
    paste(dim(U_CELL), collapse = " x "),
    "\n"
)

stopifnot(
    nrow(U_EP) == P,
    ncol(U_EP) == R_EP,
    ncol(U_CELL) == R_CELL,
    nrow(U_MODAL) == 4L,
    ncol(U_MODAL) == R_MODAL,
    K == M$K
)

if (!identical(
    rownames(U_CELL),
    M$cells
)) {
    stop("Cell order mismatch between U_CELL and manifest.")
}

if (!is.null(rownames(U_EP))) {
    expected_ep_names <- as.character(map$stage2_pair_id)

    if (!identical(
        rownames(U_EP),
        expected_ep_names
    )) {
        stop("U_EP row order does not match stage2_pair_id.")
    }
}


# ============================================================
# 2. Stage-1 aliases and offsets
# ============================================================

BA <- S$ATAC$feature_scores
BH <- S$H3K27ac$feature_scores
BR <- S$RNA$feature_scores
BD <- S$`3D`$feature_scores

VA <- S$ATAC$cell_loadings
VH <- S$H3K27ac$cell_loadings
VR <- S$RNA$cell_loadings
VD <- S$`3D`$cell_loadings

Vlist <- list(
    ATAC = VA,
    H3K27ac = VH,
    RNA = VR,
    `3D` = VD
)

offset <- list(
    ATAC = 1:K,
    H3K27ac = (K + 1):(2 * K),
    RNA = (2 * K + 1):(3 * K),
    `3D` = (3 * K + 1):(4 * K)
)

L <- 4L * K


# ============================================================
# 3. Rebuild pair-side A block
# ============================================================

build_A_block <- function(ii) {

    mp <- map[ii]

    A1 <- BA[
        mp$idx_ATAC,
        ,
        drop = FALSE
    ] * mp$w_enhancer

    A2 <- BH[
        mp$idx_H3K27ac,
        ,
        drop = FALSE
    ] * mp$w_enhancer

    A3 <- BR[
        mp$idx_RNA,
        ,
        drop = FALSE
    ] * mp$w_gene

    A4 <- BD[
        mp$idx_3D,
        ,
        drop = FALSE
    ] * mp$w_binpair

    cbind(
        A1,
        A2,
        A3,
        A4
    )
}


# ============================================================
# 4. Reverse projection of U_EP into the original 80
#    (4 modalities x 20 K) pair-side directions
#
# PAIR_PROJ[a,j] = <U_EP[,a], A_j>
# ============================================================

PAIR_PROJ <- matrix(
    0,
    nrow = R_EP,
    ncol = L
)

blocks <- split(
    seq_len(P),
    ceiling(seq_len(P) / PAIR_BLOCK)
)

cat(
    "\nComputing U_EP -> original 80 directions in ",
    length(blocks),
    " blocks...\n",
    sep = ""
)

for (ib in seq_along(blocks)) {

    ii <- blocks[[ib]]

    Ab <- build_A_block(ii)

    PAIR_PROJ <- PAIR_PROJ +
        crossprod(
            U_EP[ii, , drop = FALSE],
            Ab
        )

    if (
        ib %% 10L == 0L ||
        ib == length(blocks)
    ) {
        cat(
            "pair reverse projection ",
            ib,
            "/",
            length(blocks),
            "\n",
            sep = ""
        )
    }

    rm(Ab)
}

col_labels <- unlist(
    lapply(
        mods,
        function(m) {
            paste0(
                m,
                ":K",
                sprintf("%02d", seq_len(K))
            )
        }
    ),
    use.names = FALSE
)

rownames(PAIR_PROJ) <- sprintf(
    "EP%02d",
    seq_len(R_EP)
)

colnames(PAIR_PROJ) <- col_labels


# ============================================================
# 5. Reverse projection of U_CELL into original 80 cell-side
#    directions
#
# CELL_PROJ[[m]][b,k] = <U_CELL[,b], V_m[,k]>
#
# For tracing how U_CELL itself arose in the cell HOSVD,
# use the weighted quantity:
#
#   CELL_ORIGIN_WEIGHTED[b,m,k]
#       = ||A_mk|| * CELL_PROJ[b,m,k]
#
# because the cell HOSVD was performed on
#   [ ||A_mk|| V_mk ].
# ============================================================

CELL_PROJ <- lapply(
    Vlist,
    function(Vm) {
        crossprod(
            U_CELL,
            Vm
        )
    }
)

# Norm of every pair-side direction A_{m,k}
# This is available from PAIR_PROJ only if all R_EP directions
# were retained; they were truncated to 20. So calculate exactly
# from A in blocks.
A_NORM2 <- numeric(L)

cat("\nComputing exact norms of the original 80 pair-side directions...\n")

for (ib in seq_along(blocks)) {

    ii <- blocks[[ib]]

    Ab <- build_A_block(ii)

    A_NORM2 <- A_NORM2 +
        colSums(
            Ab^2
        )

    if (
        ib %% 10L == 0L ||
        ib == length(blocks)
    ) {
        cat(
            "A norm block ",
            ib,
            "/",
            length(blocks),
            "\n",
            sep = ""
        )
    }

    rm(Ab)
}

A_NORM <- sqrt(
    pmax(A_NORM2, 0)
)

names(A_NORM) <- col_labels

CELL_ORIGIN_WEIGHTED <- matrix(
    0,
    nrow = R_CELL,
    ncol = L
)

for (im in seq_along(mods)) {

    nm <- mods[im]
    jj <- offset[[nm]]

    CELL_ORIGIN_WEIGHTED[
        ,
        jj
    ] <- sweep(
        CELL_PROJ[[nm]],
        2,
        A_NORM[jj],
        "*"
    )
}

rownames(CELL_ORIGIN_WEIGHTED) <- sprintf(
    "CELL%02d",
    seq_len(R_CELL)
)

colnames(CELL_ORIGIN_WEIGHTED) <- col_labels


# ============================================================
# 6. Reconstruct the pre-modality-rotation core exactly
#
# CORE_PRE[a,b,k,m] =
#   PAIR_PROJ[a,m,k] * CELL_PROJ[b,m,k]
# ============================================================

CORE_PRE <- array(
    0,
    dim = c(
        R_EP,
        R_CELL,
        K,
        4L
    ),
    dimnames = list(
        EP = sprintf(
            "EP%02d",
            seq_len(R_EP)
        ),
        CELL = sprintf(
            "CELL%02d",
            seq_len(R_CELL)
        ),
        K = sprintf(
            "K%02d",
            seq_len(K)
        ),
        MODALITY = mods
    )
)

for (im in seq_along(mods)) {

    nm <- mods[im]
    jj <- offset[[nm]]

    CP <- CELL_PROJ[[nm]]

    for (k in seq_len(K)) {

        CORE_PRE[
            ,
            ,
            k,
            im
        ] <- tcrossprod(
            PAIR_PROJ[
                ,
                jj[k]
            ],
            CP[
                ,
                k
            ]
        )
    }
}

saveRDS(
    CORE_PRE,
    OUT_CORE_PRE
)


# ============================================================
# 7. Validate that CORE_PRE + U_MODAL reproduces saved CORE
# ============================================================

CORE_CHECK <- array(
    0,
    dim = dim(CORE)
)

for (r in seq_len(R_MODAL)) {

    for (im in seq_along(mods)) {

        CORE_CHECK[
            ,
            ,
            ,
            r
        ] <- CORE_CHECK[
            ,
            ,
            ,
            r
        ] +
            CORE_PRE[
                ,
                ,
                ,
                im
            ] *
            U_MODAL[
                im,
                r
            ]
    }
}

core_max_abs_error <- max(
    abs(
        CORE_CHECK -
        CORE
    )
)

core_rel_fro_error <- sqrt(
    sum(
        (
            CORE_CHECK -
            CORE
        )^2
    )
) /
    sqrt(
        sum(
            CORE^2
        )
    )

cat(
    "\nCORE reconstruction max abs error = ",
    format(core_max_abs_error, scientific = TRUE),
    "\n",
    sep = ""
)

cat(
    "CORE reconstruction relative Frobenius error = ",
    format(core_rel_fro_error, scientific = TRUE),
    "\n",
    sep = ""
)

if (core_rel_fro_error > 1e-8) {
    stop(
        "Reverse-projection reconstruction check failed."
    )
}


# ============================================================
# 8. Rank all non-K-mixed CORE[a,b,k,r] by |value|
# ============================================================

ord <- order(
    abs(as.numeric(CORE)),
    decreasing = TRUE
)

ord <- ord[
    seq_len(
        min(
            TOP_CORE,
            length(ord)
        )
    )
]

idx <- arrayInd(
    ord,
    .dim = dim(CORE)
)

top_core <- data.table(
    core_rank = seq_along(ord),
    a_EP = idx[, 1],
    b_CELL = idx[, 2],
    k = idx[, 3],
    r_MODAL = idx[, 4],
    core_value = as.numeric(CORE)[ord]
)

top_core[
    ,
    abs_core :=
        abs(core_value)
]

top_core[
    ,
    EP_component :=
        sprintf("EP%02d", a_EP)
]

top_core[
    ,
    CELL_component :=
        sprintf("CELL%02d", b_CELL)
]

top_core[
    ,
    K_component :=
        sprintf("K%02d", k)
]

top_core[
    ,
    MOD_component :=
        sprintf("MOD%02d", r_MODAL)
]

setcolorder(
    top_core,
    c(
        "core_rank",
        "a_EP",
        "b_CELL",
        "k",
        "r_MODAL",
        "EP_component",
        "CELL_component",
        "K_component",
        "MOD_component",
        "core_value",
        "abs_core"
    )
)

fwrite(
    top_core,
    OUT_TOP_CORE,
    sep = "\t"
)

cat(
    "\nTop 20 |CORE| entries:\n"
)

print(
    top_core[1:min(20L, .N)]
)


# ============================================================
# 9. Original-modality decomposition of every selected core
#
# Two different quantities are reported:
#
# formation_contribution:
#   CORE_PRE[a,b,k,m] * U_MODAL[m,r]
#   These FOUR values sum exactly to CORE[a,b,k,r].
#
# backprojected_single_core_amplitude:
#   CORE[a,b,k,r] * U_MODAL[m,r]
#   This is the amplitude of THIS ONE core element when the
#   modality HOSVD component r is mapped back to original m.
# ============================================================

modal_rows <- vector(
    "list",
    nrow(top_core)
)

for (ii in seq_len(nrow(top_core))) {

    z <- top_core[ii]

    a <- z$a_EP
    b <- z$b_CELL
    k <- z$k
    r <- z$r_MODAL
    g <- z$core_value

    pre <- CORE_PRE[
        a,
        b,
        k,

    ]

    modal_loading <- U_MODAL[
        ,
        r
    ]

    formation <- pre *
        modal_loading

    j_by_mod <- vapply(
        mods,
        function(nm) {
            offset[[nm]][k]
        },
        integer(1)
    )

    pair_projection <- vapply(
        seq_along(mods),
        function(im) {
            PAIR_PROJ[
                a,
                j_by_mod[im]
            ]
        },
        numeric(1)
    )

    cell_projection <- vapply(
        seq_along(mods),
        function(im) {
            CELL_PROJ[[mods[im]]][
                b,
                k
            ]
        },
        numeric(1)
    )

    modal_rows[[ii]] <- data.table(
        core_rank = z$core_rank,
        a_EP = a,
        b_CELL = b,
        k = k,
        r_MODAL = r,
        modality = mods,
        pair_projection = pair_projection,
        cell_projection = cell_projection,
        core_before_modality_rotation = pre,
        modality_factor_loading = modal_loading,
        formation_contribution_to_core = formation,
        backprojected_single_core_amplitude =
            g * modal_loading,
        core_value = g
    )
}

modal_decomp <- rbindlist(
    modal_rows
)

modal_check <- modal_decomp[
    ,
    .(
        reconstructed_core =
            sum(
                formation_contribution_to_core
            ),
        core_value =
            unique(
                core_value
            )
    ),
    by = core_rank
]

modal_check[
    ,
    error :=
        reconstructed_core -
        core_value
]

if (max(abs(modal_check$error)) > 1e-10) {
    stop(
        "Original modality decomposition failed."
    )
}

fwrite(
    modal_decomp,
    OUT_MODAL_DECOMP,
    sep = "\t"
)


# ============================================================
# 10. Trace selected U_EP components to original 80 directions
# ============================================================

selected_a <- sort(
    unique(
        top_core$a_EP
    )
)

ep_origin_rows <- list()
rr <- 1L

for (a in selected_a) {

    coeff <- PAIR_PROJ[
        a,

    ]

    denom <- sum(
        coeff^2
    )

    ord80 <- order(
        abs(coeff),
        decreasing = TRUE
    )

    ord80 <- ord80[
        seq_len(
            min(
                TOP_ORIGIN_80,
                L
            )
        )
    ]

    for (rank_i in seq_along(ord80)) {

        j <- ord80[rank_i]

        im <- ceiling(
            j / K
        )

        kk <- ((j - 1L) %% K) + 1L

        ep_origin_rows[[rr]] <- data.table(
            a_EP = a,
            EP_component =
                sprintf("EP%02d", a),
            origin_rank = rank_i,
            modality = mods[im],
            k = kk,
            K_component =
                sprintf("K%02d", kk),
            projection = coeff[j],
            abs_projection =
                abs(coeff[j]),
            squared_fraction =
                if (denom > 0) {
                    coeff[j]^2 / denom
                } else {
                    NA_real_
                }
        )

        rr <- rr + 1L
    }
}

ep_origin <- rbindlist(
    ep_origin_rows
)

fwrite(
    ep_origin,
    OUT_EP_ORIGIN,
    sep = "\t"
)


# ============================================================
# 11. Trace selected U_CELL components to original 80 directions
# ============================================================

selected_b <- sort(
    unique(
        top_core$b_CELL
    )
)

cell_origin_rows <- list()
rr <- 1L

for (b in selected_b) {

    coeff <- CELL_ORIGIN_WEIGHTED[
        b,

    ]

    denom <- sum(
        coeff^2
    )

    ord80 <- order(
        abs(coeff),
        decreasing = TRUE
    )

    ord80 <- ord80[
        seq_len(
            min(
                TOP_ORIGIN_80,
                L
            )
        )
    ]

    for (rank_i in seq_along(ord80)) {

        j <- ord80[rank_i]

        im <- ceiling(
            j / K
        )

        kk <- ((j - 1L) %% K) + 1L

        nm <- mods[im]

        cell_origin_rows[[rr]] <- data.table(
            b_CELL = b,
            CELL_component =
                sprintf("CELL%02d", b),
            origin_rank = rank_i,
            modality = nm,
            k = kk,
            K_component =
                sprintf("K%02d", kk),
            weighted_projection =
                coeff[j],
            abs_weighted_projection =
                abs(coeff[j]),
            squared_fraction =
                if (denom > 0) {
                    coeff[j]^2 / denom
                } else {
                    NA_real_
                },
            raw_cell_projection =
                CELL_PROJ[[nm]][
                    b,
                    kk
                ],
            pair_side_norm =
                A_NORM[j]
        )

        rr <- rr + 1L
    }
}

cell_origin <- rbindlist(
    cell_origin_rows
)

fwrite(
    cell_origin,
    OUT_CELL_ORIGIN,
    sep = "\t"
)


# ============================================================
# 12. For every selected core, return to original E-P pairs
#
# Signed score for a single core element:
#
#   core_value * U_EP[p,a]
#
# This identifies which E-P rows carry the selected EP-mode
# component, including the sign of the selected core.
#
# We output strongest positive and strongest negative E-P rows.
# ============================================================

top_ep_rows <- vector(
    "list",
    nrow(top_core)
)

for (ii in seq_len(nrow(top_core))) {

    z <- top_core[ii]

    score <- z$core_value *
        U_EP[
            ,
            z$a_EP
        ]

    pos <- order(
        score,
        decreasing = TRUE
    )[
        seq_len(
            min(
                TOP_EP_EACH_SIGN,
                length(score)
            )
        )
    ]

    neg <- order(
        score,
        decreasing = FALSE
    )[
        seq_len(
            min(
                TOP_EP_EACH_SIGN,
                length(score)
            )
        )
    ]

    idx_keep <- c(
        pos,
        neg
    )

    direction <- c(
        rep("positive", length(pos)),
        rep("negative", length(neg))
    )

    rank_within <- c(
        seq_along(pos),
        seq_along(neg)
    )

    tmp <- copy(
        map[
            idx_keep
        ]
    )

    tmp[
        ,
        core_rank :=
            z$core_rank
    ]

    tmp[
        ,
        a_EP :=
            z$a_EP
    ]

    tmp[
        ,
        b_CELL :=
            z$b_CELL
    ]

    tmp[
        ,
        k :=
            z$k
    ]

    tmp[
        ,
        r_MODAL :=
            z$r_MODAL
    ]

    tmp[
        ,
        core_value :=
            z$core_value
    ]

    tmp[
        ,
        EP_loading :=
            U_EP[
                idx_keep,
                z$a_EP
            ]
    ]

    tmp[
        ,
        signed_core_EP_score :=
            score[
                idx_keep
            ]
    ]

    tmp[
        ,
        direction :=
            direction
    ]

    tmp[
        ,
        rank_within_direction :=
            rank_within
    ]

    top_ep_rows[[ii]] <- tmp
}

top_ep <- rbindlist(
    top_ep_rows,
    fill = TRUE
)

fwrite(
    top_ep,
    OUT_TOP_EP,
    sep = "\t"
)


# ============================================================
# 13. For every selected core, return to original cells
#
# Signed score:
#
#   core_value * U_CELL[c,b]
# ============================================================

top_cell_rows <- vector(
    "list",
    nrow(top_core)
)

for (ii in seq_len(nrow(top_core))) {

    z <- top_core[ii]

    score <- z$core_value *
        U_CELL[
            ,
            z$b_CELL
        ]

    pos <- order(
        score,
        decreasing = TRUE
    )[
        seq_len(
            min(
                TOP_CELL_EACH_SIGN,
                length(score)
            )
        )
    ]

    neg <- order(
        score,
        decreasing = FALSE
    )[
        seq_len(
            min(
                TOP_CELL_EACH_SIGN,
                length(score)
            )
        )
    ]

    idx_keep <- c(
        pos,
        neg
    )

    tmp <- data.table(
        core_rank = z$core_rank,
        a_EP = z$a_EP,
        b_CELL = z$b_CELL,
        k = z$k,
        r_MODAL = z$r_MODAL,
        core_value = z$core_value,
        cell =
            rownames(U_CELL)[
                idx_keep
            ],
        CELL_loading =
            U_CELL[
                idx_keep,
                z$b_CELL
            ],
        signed_core_CELL_score =
            score[
                idx_keep
            ],
        direction = c(
            rep("positive", length(pos)),
            rep("negative", length(neg))
        ),
        rank_within_direction = c(
            seq_along(pos),
            seq_along(neg)
        )
    )

    top_cell_rows[[ii]] <- tmp
}

top_cell <- rbindlist(
    top_cell_rows
)

fwrite(
    top_cell,
    OUT_TOP_CELL,
    sep = "\t"
)


# ============================================================
# 14. Save compact RDS bundle
# ============================================================

trace_result <- list(

    top_core =
        top_core,

    original_modality_decomposition =
        modal_decomp,

    selected_EPcomponent_origin_80 =
        ep_origin,

    selected_CELLcomponent_origin_80 =
        cell_origin,

    pair_projection =
        PAIR_PROJ,

    cell_projection =
        CELL_PROJ,

    cell_origin_weighted =
        CELL_ORIGIN_WEIGHTED,

    pair_side_norm =
        A_NORM,

    core_reconstruction_max_abs_error =
        core_max_abs_error,

    core_reconstruction_relative_fro_error =
        core_rel_fro_error,

    files = list(
        top_core =
            OUT_TOP_CORE,
        original_modality_decomposition =
            OUT_MODAL_DECOMP,
        EP_origin =
            OUT_EP_ORIGIN,
        CELL_origin =
            OUT_CELL_ORIGIN,
        top_EP =
            OUT_TOP_EP,
        top_cells =
            OUT_TOP_CELL,
        core_pre =
            OUT_CORE_PRE
    )
)

saveRDS(
    trace_result,
    OUT_TRACE_RDS
)


# ============================================================
# 15. Human-readable summary of the top core element
# ============================================================

z <- top_core[1]

cat(
    "\n============================================================\n"
)

cat(
    "Top non-K-mixed CORE element\n"
)

cat(
    "============================================================\n"
)

print(z)

cat(
    "\nOriginal-modality formation decomposition:\n"
)

print(
    modal_decomp[
        core_rank == 1
    ][
        order(
            -abs(
                formation_contribution_to_core
            )
        )
    ]
)

cat(
    "\nTop original (modality,K) origins of ",
    z$EP_component,
    ":\n",
    sep = ""
)

print(
    ep_origin[
        a_EP == z$a_EP
    ]
)

cat(
    "\nTop original (modality,K) origins of ",
    z$CELL_component,
    ":\n",
    sep = ""
)

print(
    cell_origin[
        b_CELL == z$b_CELL
    ]
)

cat(
    "\nTop 10 E-P pairs for core rank 1:\n"
)

print(
    top_ep[
        core_rank == 1 &
        rank_within_direction <= 5
    ][
        ,
        .(
            direction,
            rank_within_direction,
            stage2_pair_id,
            enhancer_id,
            gene_name,
            binpair_id,
            EP_loading,
            signed_core_EP_score
        )
    ]
)

cat(
    "\nTop 10 cells for core rank 1:\n"
)

print(
    top_cell[
        core_rank == 1 &
        rank_within_direction <= 5
    ]
)


# ============================================================
# 16. Final message
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "CORE -> original modality / E-P / cell tracing completed\n"
)

cat(
    "============================================================\n"
)

cat(
    "Top CORE entries :",
    nrow(top_core),
    "\n"
)

cat(
    "CORE reconstruction relative error :",
    format(
        core_rel_fro_error,
        scientific = TRUE
    ),
    "\n"
)

cat(
    "\nSaved files:\n",
    OUT_TOP_CORE, "\n",
    OUT_MODAL_DECOMP, "\n",
    OUT_EP_ORIGIN, "\n",
    OUT_CELL_ORIGIN, "\n",
    OUT_TOP_EP, "\n",
    OUT_TOP_CELL, "\n",
    OUT_CORE_PRE, "\n",
    OUT_TRACE_RDS, "\n",
    sep = ""
)
