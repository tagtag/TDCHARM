# ============================================================
# 18_core12_CELL03_replicate_reproducibility.R
# GSE303006 / GENCODE M23
#
# Replicate-wise reproducibility analysis for the direct
# four-modality back-projection of core #12 / CELL03.
#
# Input:
#   GSE303006_M23_core12_CELL03_original4modality.rds
#
# Goal:
#   Test whether the directions observed in the pooled analysis
#   reproduce separately in R1, R2 and R3.
#
# Analyses:
#   A) Positive module:
#        RNA_logCP10K
#        ATAC_logCP10K
#        H3K27ac_logCP10K
#        3D_minuslogdistance
#        3D_distance
#
#   B) Negative module:
#        same five readouts
#
#   C) Positive - Negative module contrast:
#        RNA_z
#        ATAC_z
#        H3K27ac_z
#        3D_z
#
# For each replicate:
#   - n positive / negative CELL03 cells
#   - mean in each group
#   - delta mean = positive - negative
#   - median difference
#   - rank-biserial effect size
#   - Wilcoxon p value
#
# Across R1/R2/R3:
#   - sign concordance
#   - number of replicates with same direction
#   - median replicate effect
#   - fixed-direction Stouffer combined p (sign-aware)
#
# NOTE:
#   This is within-dataset reproducibility across biological
#   replicates, not an independent external validation.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

needed <- c(
    "data.table"
)

missing <- needed[
    !vapply(
        needed,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]

if (length(missing) > 0L) {
    stop(
        "Missing package(s): ",
        paste(missing, collapse = ", ")
    )
}

library(data.table)


# ============================================================
# 1. Files / settings
# ============================================================

INPUT_RDS <-
    "GSE303006_M23_core12_CELL03_original4modality.rds"

OUT_PREFIX <-
    "GSE303006_M23_core12_CELL03_replicate_reproducibility"

REPLICATES <- c(
    "R1",
    "R2",
    "R3"
)


# ============================================================
# 2. Input
# ============================================================

if (!file.exists(INPUT_RDS)) {
    stop(
        "Missing input file: ",
        INPUT_RDS
    )
}

obj <- readRDS(
    INPUT_RDS
)

if (!"module_cell_scores" %in% names(obj)) {
    stop(
        "Input RDS lacks module_cell_scores."
    )
}

dt <- as.data.table(
    copy(
        obj$module_cell_scores
    )
)

required_cols <- c(
    "cell",
    "CELL03_loading",
    "group",
    "replicate",

    "POS_RNA_logCP10K",
    "POS_ATAC_logCP10K",
    "POS_H3K27ac_logCP10K",
    "POS_3D_minuslogdistance",
    "POS_3D_distance",

    "NEG_RNA_logCP10K",
    "NEG_ATAC_logCP10K",
    "NEG_H3K27ac_logCP10K",
    "NEG_3D_minuslogdistance",
    "NEG_3D_distance",

    "contrast_RNA_z",
    "contrast_ATAC_z",
    "contrast_H3K27ac_z",
    "contrast_3D_z"
)

missing_cols <- setdiff(
    required_cols,
    names(dt)
)

if (length(missing_cols) > 0L) {
    stop(
        "Missing columns in module_cell_scores: ",
        paste(
            missing_cols,
            collapse = ", "
        )
    )
}

dt[
    ,
    group :=
        factor(
            as.character(group),
            levels = c(
                "negative",
                "positive"
            )
        )
]

dt[
    ,
    replicate :=
        as.character(
            replicate
        )
]

cat(
    "Selected cells total :",
    nrow(dt),
    "\n"
)

cat(
    "\nCells by replicate and CELL03 sign:\n"
)

print(
    table(
        dt$replicate,
        dt$group,
        useNA = "ifany"
    )
)


# ============================================================
# 3. Helper functions
# ============================================================

rank_biserial <- function(
    x,
    y
) {

    x <- x[
        is.finite(x)
    ]

    y <- y[
        is.finite(y)
    ]

    n1 <- length(x)
    n2 <- length(y)

    if (
        n1 == 0L ||
        n2 == 0L
    ) {
        return(
            NA_real_
        )
    }

    rr <- rank(
        c(
            x,
            y
        ),
        ties.method = "average"
    )

    U <- sum(
        rr[
            seq_len(n1)
        ]
    ) -
        n1 *
        (
            n1 + 1
        ) /
        2

    2 *
        U /
        (
            n1 *
            n2
        ) -
        1
}


one_replicate_comparison <- function(
    x,
    group
) {

    ok <-
        is.finite(x) &
        !is.na(group)

    x <- x[ok]
    group <- group[ok]

    xp <- x[
        group == "positive"
    ]

    xn <- x[
        group == "negative"
    ]

    if (
        length(xp) == 0L ||
        length(xn) == 0L
    ) {

        return(
            data.table(
                n_positive =
                    length(xp),

                n_negative =
                    length(xn),

                mean_positive =
                    NA_real_,

                mean_negative =
                    NA_real_,

                delta_mean_pos_minus_neg =
                    NA_real_,

                median_positive =
                    NA_real_,

                median_negative =
                    NA_real_,

                delta_median_pos_minus_neg =
                    NA_real_,

                rank_biserial =
                    NA_real_,

                wilcox_p =
                    NA_real_
            )
        )
    }

    wt <- suppressWarnings(
        wilcox.test(
            xp,
            xn,
            exact = FALSE
        )
    )

    data.table(
        n_positive =
            length(xp),

        n_negative =
            length(xn),

        mean_positive =
            mean(
                xp
            ),

        mean_negative =
            mean(
                xn
            ),

        delta_mean_pos_minus_neg =
            mean(xp) -
            mean(xn),

        median_positive =
            median(
                xp
            ),

        median_negative =
            median(
                xn
            ),

        delta_median_pos_minus_neg =
            median(xp) -
            median(xn),

        rank_biserial =
            rank_biserial(
                xp,
                xn
            ),

        wilcox_p =
            wt$p.value
    )
}


signed_stouffer <- function(
    p,
    effect,
    weights = NULL
) {

    ok <-
        is.finite(p) &
        p > 0 &
        p <= 1 &
        is.finite(effect) &
        effect != 0

    p <- p[ok]
    effect <- effect[ok]

    if (is.null(weights)) {
        weights <- rep(
            1,
            length(p)
        )
    } else {
        weights <- weights[ok]
    }

    if (length(p) == 0L) {
        return(
            c(
                Z = NA_real_,
                p_two_sided = NA_real_
            )
        )
    }

    # Two-sided p -> signed z.
    z_abs <- qnorm(
        1 -
        p /
        2
    )

    z <- sign(effect) *
        z_abs

    Z <- sum(
        weights *
        z
    ) /
        sqrt(
            sum(
                weights^2
            )
        )

    p_comb <- 2 *
        pnorm(
            -abs(Z)
        )

    c(
        Z = Z,
        p_two_sided = p_comb
    )
}


direction_label <- function(
    x
) {

    fifelse(
        x > 0,
        "positive>negative",
        fifelse(
            x < 0,
            "positive<negative",
            "equal"
        )
    )
}


# ============================================================
# 4. Define analyses
# ============================================================

module_specs <- rbind(

    data.table(
        analysis_type =
            "positive_module",

        modality =
            c(
                "RNA_logCP10K",
                "ATAC_logCP10K",
                "H3K27ac_logCP10K",
                "3D_minuslogdistance",
                "3D_distance"
            ),

        column =
            c(
                "POS_RNA_logCP10K",
                "POS_ATAC_logCP10K",
                "POS_H3K27ac_logCP10K",
                "POS_3D_minuslogdistance",
                "POS_3D_distance"
            )
    ),

    data.table(
        analysis_type =
            "negative_module",

        modality =
            c(
                "RNA_logCP10K",
                "ATAC_logCP10K",
                "H3K27ac_logCP10K",
                "3D_minuslogdistance",
                "3D_distance"
            ),

        column =
            c(
                "NEG_RNA_logCP10K",
                "NEG_ATAC_logCP10K",
                "NEG_H3K27ac_logCP10K",
                "NEG_3D_minuslogdistance",
                "NEG_3D_distance"
            )
    ),

    data.table(
        analysis_type =
            "POS_minus_NEG_contrast",

        modality =
            c(
                "RNA_z",
                "ATAC_z",
                "H3K27ac_z",
                "3D_z"
            ),

        column =
            c(
                "contrast_RNA_z",
                "contrast_ATAC_z",
                "contrast_H3K27ac_z",
                "contrast_3D_z"
            )
    )
)


# ============================================================
# 5. Replicate-wise comparisons
# ============================================================

replicate_rows <- list()

rr <- 1L

for (i in seq_len(
    nrow(
        module_specs
    )
)) {

    sp <- module_specs[i]

    for (rep in REPLICATES) {

        drep <- dt[
            replicate == rep
        ]

        cmp <- one_replicate_comparison(
            drep[[
                    sp$column
                ]
            ],
            drep$group
        )

        cmp[
            ,
            analysis_type :=
                sp$analysis_type
        ]

        cmp[
            ,
            modality :=
                sp$modality
        ]

        cmp[
            ,
            score_column :=
                sp$column
        ]

        cmp[
            ,
            replicate :=
                rep
        ]

        cmp[
            ,
            direction_delta_mean :=
                direction_label(
                    delta_mean_pos_minus_neg
                )
        ]

        cmp[
            ,
            direction_rank_biserial :=
                direction_label(
                    rank_biserial
                )
        ]

        replicate_rows[[rr]] <-
            cmp

        rr <- rr + 1L
    }
}

replicate_results <- rbindlist(
    replicate_rows,
    fill = TRUE
)

replicate_results[
    ,
    wilcox_FDR_all_tests :=
        p.adjust(
            wilcox_p,
            method = "BH"
        )
]

replicate_results[
    ,
    wilcox_FDR_within_analysis :=
        p.adjust(
            wilcox_p,
            method = "BH"
        ),
    by = analysis_type
]

setcolorder(
    replicate_results,
    c(
        "analysis_type",
        "modality",
        "replicate",
        "n_positive",
        "n_negative",
        "mean_positive",
        "mean_negative",
        "delta_mean_pos_minus_neg",
        "rank_biserial",
        "direction_delta_mean",
        "wilcox_p",
        "wilcox_FDR_within_analysis",
        "wilcox_FDR_all_tests",
        setdiff(
            names(
                replicate_results
            ),
            c(
                "analysis_type",
                "modality",
                "replicate",
                "n_positive",
                "n_negative",
                "mean_positive",
                "mean_negative",
                "delta_mean_pos_minus_neg",
                "rank_biserial",
                "direction_delta_mean",
                "wilcox_p",
                "wilcox_FDR_within_analysis",
                "wilcox_FDR_all_tests"
            )
        )
    )
)

fwrite(
    replicate_results,
    paste0(
        OUT_PREFIX,
        "_replicate_effects.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 6. Reproducibility summary across R1/R2/R3
# ============================================================

repro_summary <- replicate_results[
    ,
    {

        delta <- delta_mean_pos_minus_neg
        rb <- rank_biserial
        p <- wilcox_p

        valid <- is.finite(
            delta
        )

        delta_valid <- delta[
            valid
        ]

        rb_valid <- rb[
            valid
        ]

        p_valid <- p[
            valid
        ]

        npos <- n_positive[
            valid
        ]

        nneg <- n_negative[
            valid
        ]

        weights <- sqrt(
            pmax(
                npos +
                nneg,
                1
            )
        )

        st <- signed_stouffer(
            p =
                p_valid,

            effect =
                rb_valid,

            weights =
                weights
        )

        signs <- sign(
            delta_valid
        )

        nonzero_signs <- signs[
            signs != 0
        ]

        all_same_sign <-
            length(nonzero_signs) > 0L &&
            length(
                unique(
                    nonzero_signs
                )
            ) == 1L

        majority_sign <- if (
            length(nonzero_signs) == 0L
        ) {

            0L

        } else {

            s <- sum(
                nonzero_signs
            )

            if (s > 0) {
                1L
            } else if (s < 0) {
                -1L
            } else {
                0L
            }
        }

        n_same_majority <- if (
            majority_sign == 0L
        ) {

            0L

        } else {

            sum(
                signs ==
                    majority_sign,
                na.rm = TRUE
            )
        }

        list(
            n_replicates_available =
                sum(valid),

            n_same_direction_majority =
                n_same_majority,

            all_three_same_direction =
                sum(valid) == 3L &&
                all_same_sign,

            majority_direction =
                if (
                    majority_sign > 0
                ) {
                    "positive>negative"
                } else if (
                    majority_sign < 0
                ) {
                    "positive<negative"
                } else {
                    "mixed"
                },

            median_delta_mean =
                median(
                    delta_valid,
                    na.rm = TRUE
                ),

            min_delta_mean =
                min(
                    delta_valid,
                    na.rm = TRUE
                ),

            max_delta_mean =
                max(
                    delta_valid,
                    na.rm = TRUE
                ),

            median_rank_biserial =
                median(
                    rb_valid,
                    na.rm = TRUE
                ),

            min_rank_biserial =
                min(
                    rb_valid,
                    na.rm = TRUE
                ),

            max_rank_biserial =
                max(
                    rb_valid,
                    na.rm = TRUE
                ),

            signed_Stouffer_Z =
                unname(
                    st[
                        "Z"
                    ]
                ),

            signed_Stouffer_p =
                unname(
                    st[
                        "p_two_sided"
                    ]
                )
        )
    },
    by = .(
        analysis_type,
        modality,
        score_column
    )
]

repro_summary[
    ,
    signed_Stouffer_FDR :=
        p.adjust(
            signed_Stouffer_p,
            method = "BH"
        )
]

setorder(
    repro_summary,
    analysis_type,
    modality
)

fwrite(
    repro_summary,
    paste0(
        OUT_PREFIX,
        "_summary.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 7. Direction matrix for easy inspection
# ============================================================

direction_matrix <- dcast(
    replicate_results,
    analysis_type + modality ~ replicate,
    value.var = "delta_mean_pos_minus_neg"
)

fwrite(
    direction_matrix,
    paste0(
        OUT_PREFIX,
        "_delta_mean_matrix.tsv"
    ),
    sep = "\t"
)

rb_matrix <- dcast(
    replicate_results,
    analysis_type + modality ~ replicate,
    value.var = "rank_biserial"
)

fwrite(
    rb_matrix,
    paste0(
        OUT_PREFIX,
        "_rank_biserial_matrix.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 8. Plot replicate effect sizes
# ============================================================

pdf(
    paste0(
        OUT_PREFIX,
        "_rank_biserial_by_replicate.pdf"
    ),
    width = 10,
    height = 7
)

analysis_levels <- unique(
    replicate_results$analysis_type
)

for (aa in analysis_levels) {

    dd <- replicate_results[
        analysis_type == aa
    ]

    modalities <- unique(
        dd$modality
    )

    mat <- matrix(
        NA_real_,
        nrow = length(modalities),
        ncol = length(REPLICATES),
        dimnames = list(
            modalities,
            REPLICATES
        )
    )

    for (m in modalities) {

        zz <- dd[
            modality == m
        ]

        mat[
            m,
            zz$replicate
        ] <- zz$rank_biserial
    }

    lim <- max(
        abs(
            mat
        ),
        na.rm = TRUE
    )

    if (!is.finite(lim) ||
        lim == 0) {
        lim <- 1
    }

    matplot(
        t(
            mat
        ),
        type = "b",
        pch = seq_len(
            nrow(mat)
        ),
        lty = seq_len(
            nrow(mat)
        ),
        xaxt = "n",
        xlab = "Replicate",
        ylab = "Rank-biserial effect (CELL03 positive - negative)",
        ylim = c(
            -lim,
            lim
        ),
        main = aa
    )

    axis(
        1,
        at = seq_along(
            REPLICATES
        ),
        labels = REPLICATES
    )

    abline(
        h = 0,
        lty = 2
    )

    legend(
        "topright",
        legend = rownames(mat),
        pch = seq_len(
            nrow(mat)
        ),
        lty = seq_len(
            nrow(mat)
        ),
        cex = 0.8,
        bty = "n"
    )
}

dev.off()


# ============================================================
# 9. Heatmap-like signed effect plot
# ============================================================

effect_labels <- paste(
    replicate_results$analysis_type,
    replicate_results$modality,
    sep = " | "
)

effect_levels <- unique(
    effect_labels
)

effect_mat <- matrix(
    NA_real_,
    nrow = length(effect_levels),
    ncol = length(REPLICATES),
    dimnames = list(
        effect_levels,
        REPLICATES
    )
)

for (i in seq_len(
    nrow(
        replicate_results
    )
)) {

    lab <- paste(
        replicate_results$analysis_type[i],
        replicate_results$modality[i],
        sep = " | "
    )

    effect_mat[
        lab,
        replicate_results$replicate[i]
    ] <- replicate_results$rank_biserial[i]
}

clip <- 1

effect_mat_plot <- pmax(
    pmin(
        effect_mat,
        clip
    ),
    -clip
)

cols <- hcl.colors(
    101,
    "Blue-Red 3",
    rev = TRUE
)

pdf(
    paste0(
        OUT_PREFIX,
        "_effect_heatmap.pdf"
    ),
    width = 8,
    height = max(
        8,
        0.35 *
            nrow(effect_mat_plot) +
            2
    )
)

par(
    mar = c(
        5,
        18,
        4,
        2
    )
)

image(
    x = seq_len(
        ncol(
            effect_mat_plot
        )
    ),
    y = seq_len(
        nrow(
            effect_mat_plot
        )
    ),
    z = t(
        effect_mat_plot[
            nrow(
                effect_mat_plot
            ):1,
            ,
            drop = FALSE
        ]
    ),
    col = cols,
    zlim = c(
        -1,
        1
    ),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Replicate-wise rank-biserial effects"
)

axis(
    1,
    at = seq_len(
        ncol(
            effect_mat_plot
        )
    ),
    labels = colnames(
        effect_mat_plot
    )
)

axis(
    2,
    at = seq_len(
        nrow(
            effect_mat_plot
        )
    ),
    labels = rev(
        rownames(
            effect_mat_plot
        )
    ),
    las = 2,
    cex.axis = 0.65
)

box()

dev.off()


# ============================================================
# 10. Console summary
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Replicate-wise reproducibility completed\n"
)

cat(
    "============================================================\n"
)

cat(
    "\nReplicate effects:\n"
)

print(
    replicate_results[
        ,
        .(
            analysis_type,
            modality,
            replicate,
            n_positive,
            n_negative,
            delta_mean_pos_minus_neg,
            rank_biserial,
            wilcox_p,
            wilcox_FDR_within_analysis
        )
    ]
)

cat(
    "\nReproducibility summary across R1/R2/R3:\n"
)

print(
    repro_summary[
        ,
        .(
            analysis_type,
            modality,
            n_replicates_available,
            n_same_direction_majority,
            all_three_same_direction,
            majority_direction,
            median_delta_mean,
            median_rank_biserial,
            min_rank_biserial,
            max_rank_biserial,
            signed_Stouffer_Z,
            signed_Stouffer_p,
            signed_Stouffer_FDR
        )
    ]
)


# ============================================================
# 11. Highlight key four-modality contrast
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "POS - NEG module contrast: key reproducibility table\n"
)

cat(
    "============================================================\n"
)

print(
    replicate_results[
        analysis_type ==
            "POS_minus_NEG_contrast",
        .(
            modality,
            replicate,
            n_positive,
            n_negative,
            delta_mean_pos_minus_neg,
            rank_biserial,
            wilcox_p
        )
    ]
)

cat(
    "\nContrast sign concordance:\n"
)

print(
    repro_summary[
        analysis_type ==
            "POS_minus_NEG_contrast",
        .(
            modality,
            all_three_same_direction,
            majority_direction,
            median_rank_biserial,
            signed_Stouffer_p,
            signed_Stouffer_FDR
        )
    ]
)


# ============================================================
# 12. Save RDS
# ============================================================

result <- list(

    source =
        INPUT_RDS,

    selected_cells =
        dt,

    replicate_effects =
        replicate_results,

    reproducibility_summary =
        repro_summary,

    delta_mean_matrix =
        direction_matrix,

    rank_biserial_matrix =
        rb_matrix
)

saveRDS(
    result,
    paste0(
        OUT_PREFIX,
        ".rds"
    )
)


# ============================================================
# 13. Final file list
# ============================================================

cat(
    "\nSaved files:\n",
    paste0(OUT_PREFIX, "_replicate_effects.tsv\n"),
    paste0(OUT_PREFIX, "_summary.tsv\n"),
    paste0(OUT_PREFIX, "_delta_mean_matrix.tsv\n"),
    paste0(OUT_PREFIX, "_rank_biserial_matrix.tsv\n"),
    paste0(OUT_PREFIX, "_rank_biserial_by_replicate.pdf\n"),
    paste0(OUT_PREFIX, "_effect_heatmap.pdf\n"),
    paste0(OUT_PREFIX, ".rds\n"),
    sep = ""
)
