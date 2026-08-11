# ============================================================
# 19_top100_CELLcomponent_celltype_vs_replicate_effect.R
# GSE303006 / GENCODE M23
#
# Generalization figure:
# Are CELL components represented among the top 100 stage-2
# core elements more strongly associated with cell type than
# with biological replicate?
#
# Input:
#   GSE303006_M23_stage2_top100_core_cell_annotation.tsv
#
# One point = one DISTINCT CELL component appearing in top 100.
#
# x-axis:
#   replicate Kruskal-Wallis epsilon^2
#
# y-axis:
#   cell-type Kruskal-Wallis epsilon^2
#
# Point size:
#   number of top-100 core elements containing that CELL component
#
# The diagonal y=x separates:
#   above line -> cell-type effect > replicate effect
#   below line -> replicate effect > cell-type effect
#
# Outputs:
#   - PDF and PNG publication-style figure
#   - component-level summary table
#   - concise statistics table
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 0. Files
# ------------------------------------------------------------

INPUT_FILE <-
    "GSE303006_M23_stage2_top100_core_cell_annotation.tsv"

OUT_COMPONENTS <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate.tsv"

OUT_STATS <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate_stats.tsv"

OUT_PDF <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate.pdf"

OUT_PNG <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate.png"

OUT_RDS <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate.rds"


# ------------------------------------------------------------
# 1. Read top-100 core annotation
# ------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
    stop(
        "Missing input file: ",
        INPUT_FILE
    )
}

x <- fread(
    INPUT_FILE
)

required <- c(
    "core_rank",
    "b_CELL",
    "CELL_component",
    "abs_core",
    "celltype_KW_epsilon2",
    "celltype_KW_FDR",
    "replicate_KW_epsilon2",
    "replicate_KW_FDR",
    "positive_top50_celltype",
    "negative_top50_celltype"
)

missing <- setdiff(
    required,
    names(x)
)

if (length(missing) > 0L) {
    stop(
        "Missing required columns: ",
        paste(
            missing,
            collapse = ", "
        )
    )
}

if (nrow(x) != 100L) {
    warning(
        "Expected 100 core rows, found ",
        nrow(x)
    )
}


# ------------------------------------------------------------
# 2. Verify effect sizes are constant within CELL component
# ------------------------------------------------------------

check <- x[
    ,
    .(
        n_celltype_epsilon =
            uniqueN(
                celltype_KW_epsilon2
            ),

        n_replicate_epsilon =
            uniqueN(
                replicate_KW_epsilon2
            )
    ),
    by = .(
        b_CELL,
        CELL_component
    )
]

if (
    any(
        check$n_celltype_epsilon != 1L
    ) ||
    any(
        check$n_replicate_epsilon != 1L
    )
) {
    stop(
        "Effect sizes are not constant within CELL component."
    )
}


# ------------------------------------------------------------
# 3. Collapse top 100 cores to DISTINCT CELL components
# ------------------------------------------------------------

comp <- x[
    ,
    .(
        n_top100_cores =
            .N,

        best_core_rank =
            min(
                core_rank
            ),

        max_abs_core =
            max(
                abs_core
            ),

        celltype_epsilon2 =
            first(
                celltype_KW_epsilon2
            ),

        replicate_epsilon2 =
            first(
                replicate_KW_epsilon2
            ),

        celltype_FDR =
            first(
                celltype_KW_FDR
            ),

        replicate_FDR =
            first(
                replicate_KW_FDR
            ),

        positive_celltype =
            first(
                positive_top50_celltype
            ),

        negative_celltype =
            first(
                negative_top50_celltype
            )
    ),
    by = .(
        b_CELL,
        CELL_component
    )
]

comp[
    ,
    epsilon2_difference :=
        celltype_epsilon2 -
        replicate_epsilon2
]

comp[
    ,
    celltype_gt_replicate :=
        celltype_epsilon2 >
        replicate_epsilon2
]

comp[
    ,
    effect_ratio :=
        (celltype_epsilon2 + 1e-12) /
        (replicate_epsilon2 + 1e-12)
]

setorder(
    comp,
    -celltype_epsilon2,
    replicate_epsilon2
)

fwrite(
    comp,
    OUT_COMPONENTS,
    sep = "\t"
)


# ------------------------------------------------------------
# 4. Overall descriptive statistics
# ------------------------------------------------------------

n_components <- nrow(
    comp
)

n_celltype_gt <- sum(
    comp$celltype_gt_replicate
)

fraction_celltype_gt <-
    n_celltype_gt /
    n_components

# Top-100 core entries whose CELL component is cell-type dominated.
celltype_dominated_b <- comp[
    celltype_gt_replicate == TRUE,
    b_CELL
]

n_core_entries_celltype_gt <- x[
    b_CELL %in%
        celltype_dominated_b,
    .N
]

fraction_core_entries_celltype_gt <-
    n_core_entries_celltype_gt /
    nrow(x)

# One-sided sign test:
# Under a null with no preferred direction, P(celltype > replicate)=0.5.
sign_test <- binom.test(
    x =
        n_celltype_gt,
    n =
        n_components,
    p =
        0.5,
    alternative =
        "greater"
)

# Paired Wilcoxon across DISTINCT CELL components.
# This is supplementary/descriptive because HOSVD components
# are not conventional independent biological observations.
wilcox_test <- suppressWarnings(
    wilcox.test(
        comp$celltype_epsilon2,
        comp$replicate_epsilon2,
        paired = TRUE,
        alternative = "greater",
        exact = FALSE
    )
)

stats <- data.table(
    metric = c(
        "n_distinct_CELL_components",
        "n_components_celltype_gt_replicate",
        "fraction_components_celltype_gt_replicate",
        "n_top100_core_entries_celltype_gt_replicate",
        "fraction_top100_core_entries_celltype_gt_replicate",
        "median_celltype_epsilon2",
        "median_replicate_epsilon2",
        "median_epsilon2_difference",
        "min_epsilon2_difference",
        "max_epsilon2_difference",
        "one_sided_sign_test_p",
        "paired_Wilcoxon_greater_p"
    ),

    value = c(
        n_components,
        n_celltype_gt,
        fraction_celltype_gt,
        n_core_entries_celltype_gt,
        fraction_core_entries_celltype_gt,
        median(
            comp$celltype_epsilon2
        ),
        median(
            comp$replicate_epsilon2
        ),
        median(
            comp$epsilon2_difference
        ),
        min(
            comp$epsilon2_difference
        ),
        max(
            comp$epsilon2_difference
        ),
        sign_test$p.value,
        wilcox_test$p.value
    )
)

fwrite(
    stats,
    OUT_STATS,
    sep = "\t"
)


# ------------------------------------------------------------
# 5. Figure helper
# ------------------------------------------------------------

draw_figure <- function() {

    # Point size encodes frequency among top 100 cores.
    # sqrt scaling keeps large-frequency components manageable.
    cex_point <-
        1.4 +
        0.45 *
        sqrt(
            comp$n_top100_cores
        )

    # Plot limits include all points and enough space for labels.
    lim_max <- max(
        c(
            comp$celltype_epsilon2,
            comp$replicate_epsilon2
        ),
        na.rm = TRUE
    )

    lim_max <- max(
        0.85,
        lim_max * 1.05
    )

    par(
        mar = c(
            5.2,
            5.5,
            3.3,
            2.0
        ),
        mgp = c(
            3.2,
            0.9,
            0
        )
    )

    plot(
        comp$replicate_epsilon2,
        comp$celltype_epsilon2,
        xlim = c(
            0,
            lim_max
        ),
        ylim = c(
            0,
            lim_max
        ),
        xlab =
            expression(
                "Replicate effect  " *
                epsilon^2
            ),
        ylab =
            expression(
                "Cell-type effect  " *
                epsilon^2
            ),
        pch = 21,
        bg = "white",
        col = "black",
        lwd = 1.5,
        cex = cex_point,
        asp = 1,
        main =
            "Top-100 core CELL components are dominated by cell-type effects"
    )

    # y = x reference.
    abline(
        a = 0,
        b = 1,
        lty = 2,
        lwd = 1.5,
        col = "grey45"
    )

    # Light guide line near the observed replicate-effect range.
    grid(
        nx = NULL,
        ny = NULL,
        lty = 3,
        col = "grey90"
    )

    # Redraw points above grid.
    points(
        comp$replicate_epsilon2,
        comp$celltype_epsilon2,
        pch = 21,
        bg = "white",
        col = "black",
        lwd = 1.5,
        cex = cex_point
    )

    # Labels. Include number of top-100 cores in parentheses.
    labels <- paste0(
        comp$CELL_component,
        " (",
        comp$n_top100_cores,
        ")"
    )

    # Because x values are small and points are well separated mostly by y,
    # put labels just to the right of the point.
    text(
        comp$replicate_epsilon2,
        comp$celltype_epsilon2,
        labels = labels,
        pos = 4,
        offset = 0.55,
        cex = 0.75
    )

    # Add compact summary directly inside figure.
    summary_text <- paste0(
        n_celltype_gt,
        "/",
        n_components,
        " distinct CELL components above y=x\n",
        n_core_entries_celltype_gt,
        "/",
        nrow(x),
        " top-core entries use these components\n",
        "Sign test P = ",
        format(
            sign_test$p.value,
            scientific = TRUE,
            digits = 3
        )
    )

    legend(
        "bottomright",
        legend = summary_text,
        bty = "n",
        cex = 0.82,
        text.col = "black"
    )

    # Point-size key.
    size_examples <- sort(
        unique(
            c(
                min(
                    comp$n_top100_cores
                ),
                round(
                    median(
                        comp$n_top100_cores
                    )
                ),
                max(
                    comp$n_top100_cores
                )
            )
        )
    )

    size_cex <-
        1.4 +
        0.45 *
        sqrt(
            size_examples
        )

    legend(
        "topright",
        legend = paste0(
            size_examples,
            " top-100 cores"
        ),
        pt.cex = size_cex,
        pch = 21,
        pt.bg = "white",
        col = "black",
        title = "Point size",
        bty = "n",
        cex = 0.72
    )
}


# ------------------------------------------------------------
# 6. Save PDF
# ------------------------------------------------------------

pdf(
    OUT_PDF,
    width = 7.5,
    height = 7.0,
    useDingbats = FALSE
)

draw_figure()

dev.off()


# ------------------------------------------------------------
# 7. Save PNG
# ------------------------------------------------------------

png(
    OUT_PNG,
    width = 1800,
    height = 1680,
    res = 240
)

draw_figure()

dev.off()


# ------------------------------------------------------------
# 8. Console report
# ------------------------------------------------------------

cat(
    "\n============================================================\n"
)

cat(
    "Top-100 CELL component generalization analysis\n"
)

cat(
    "============================================================\n"
)

cat(
    "Distinct CELL components represented in top 100 cores : ",
    n_components,
    "\n",
    sep = ""
)

cat(
    "CELL components with cell-type epsilon^2 > replicate epsilon^2 : ",
    n_celltype_gt,
    " / ",
    n_components,
    " (",
    sprintf(
        "%.1f%%",
        100 *
            fraction_celltype_gt
    ),
    ")\n",
    sep = ""
)

cat(
    "Top-100 core entries using cell-type-dominated components : ",
    n_core_entries_celltype_gt,
    " / ",
    nrow(x),
    " (",
    sprintf(
        "%.1f%%",
        100 *
            fraction_core_entries_celltype_gt
    ),
    ")\n",
    sep = ""
)

cat(
    "\nMedian cell-type epsilon^2 : ",
    median(
        comp$celltype_epsilon2
    ),
    "\n",
    sep = ""
)

cat(
    "Median replicate epsilon^2 : ",
    median(
        comp$replicate_epsilon2
    ),
    "\n",
    sep = ""
)

cat(
    "Median difference          : ",
    median(
        comp$epsilon2_difference
    ),
    "\n",
    sep = ""
)

cat(
    "\nOne-sided sign test P      : ",
    sign_test$p.value,
    "\n",
    sep = ""
)

cat(
    "Paired Wilcoxon P          : ",
    wilcox_test$p.value,
    "\n",
    sep = ""
)

cat(
    "\nComponent-level table:\n"
)

print(
    comp[
        ,
        .(
            CELL_component,
            n_top100_cores,
            best_core_rank,
            positive_celltype,
            negative_celltype,
            celltype_epsilon2,
            replicate_epsilon2,
            epsilon2_difference,
            effect_ratio
        )
    ]
)

cat(
    "\nSaved files:\n",
    OUT_COMPONENTS, "\n",
    OUT_STATS, "\n",
    OUT_PDF, "\n",
    OUT_PNG, "\n",
    OUT_RDS, "\n",
    sep = ""
)


# ------------------------------------------------------------
# 9. Save RDS
# ------------------------------------------------------------

saveRDS(
    list(
        top100_core_rows =
            x,

        component_summary =
            comp,

        statistics =
            stats,

        sign_test =
            sign_test,

        paired_wilcoxon =
            wilcox_test
    ),
    OUT_RDS
)
