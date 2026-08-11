# ============================================================
# 19b_top100_CELLcomponent_publication_figure.R
# Publication-quality revision of:
#   cell-type effect epsilon^2 vs replicate effect epsilon^2
#
# Input:
#   GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate.tsv
#
# Main changes from the first figure:
#   - manual label positions to avoid overlap
#   - CELL03 highlighted as the detailed case study
#   - other components shown neutrally
#   - point-size legend removed
#   - compact summary placed safely inside empty lower-right region
#   - regions above/below y=x explicitly labelled
#   - both PDF and high-resolution PNG output
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 0. Files
# ------------------------------------------------------------

INPUT_FILE <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate.tsv"

OUT_PDF <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate_publication.pdf"

OUT_PNG <-
    "GSE303006_M23_top100_CELLcomponents_celltype_vs_replicate_publication.png"


# ------------------------------------------------------------
# 1. Read data
# ------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
    stop(
        "Missing input file: ",
        INPUT_FILE
    )
}

comp <- fread(
    INPUT_FILE
)

required <- c(
    "CELL_component",
    "n_top100_cores",
    "celltype_epsilon2",
    "replicate_epsilon2"
)

missing <- setdiff(
    required,
    names(comp)
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


# ------------------------------------------------------------
# 2. Summary statistics
# ------------------------------------------------------------

n_components <- nrow(
    comp
)

n_above <- sum(
    comp$celltype_epsilon2 >
        comp$replicate_epsilon2
)

sign_test <- binom.test(
    n_above,
    n_components,
    p = 0.5,
    alternative = "greater"
)


# ------------------------------------------------------------
# 3. Manual label offsets
#
# Offsets are in data-coordinate units.
# These are tuned for the observed 12 CELL components.
#
# If a future analysis changes component positions substantially,
# only this table needs adjustment.
# ------------------------------------------------------------

label_offsets <- data.table(

    CELL_component = c(
        "CELL03",
        "CELL07",
        "CELL05",
        "CELL06",
        "CELL04",
        "CELL08",
        "CELL09",
        "CELL02",
        "CELL10",
        "CELL11",
        "CELL01",
        "CELL13"
    ),

    dx = c(
         0.020,
         0.020,
         0.025,
         0.025,
         0.020,
         0.020,
         0.018,
         0.020,
         0.020,
         0.020,
         0.020,
         0.020
    ),

    dy = c(
         0.000,
         0.015,
        -0.010,
        -0.020,
         0.018,
        -0.015,
         0.000,
         0.000,
         0.012,
        -0.012,
         0.012,
        -0.012
    )
)

comp <- merge(
    comp,
    label_offsets,
    by = "CELL_component",
    all.x = TRUE,
    sort = FALSE
)

comp[
    is.na(dx),
    dx := 0.02
]

comp[
    is.na(dy),
    dy := 0
]

comp[
    ,
    label_x :=
        replicate_epsilon2 +
        dx
]

comp[
    ,
    label_y :=
        celltype_epsilon2 +
        dy
]


# ------------------------------------------------------------
# 4. Point size
#
# Area approximately scales with frequency.
# ------------------------------------------------------------

point_cex <-
    1.25 +
    0.32 *
    sqrt(
        comp$n_top100_cores
    )


# ------------------------------------------------------------
# 5. Figure drawing function
# ------------------------------------------------------------

draw_publication_figure <- function() {

    xlim <- c(
        -0.015,
        0.86
    )

    ylim <- c(
        -0.015,
        0.87
    )

    par(
        mar = c(
            5.4,
            5.8,
            2.3,
            1.5
        ),
        mgp = c(
            3.2,
            0.9,
            0
        ),
        las = 1
    )

    # Empty plot first.
    plot(
        NA,
        xlim = xlim,
        ylim = ylim,
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
        xaxt = "n",
        yaxt = "n",
        bty = "l"
    )

    axis(
        1,
        at = seq(
            0,
            0.8,
            by = 0.2
        ),
        labels = sprintf(
            "%.1f",
            seq(
                0,
                0.8,
                by = 0.2
            )
        )
    )

    axis(
        2,
        at = seq(
            0,
            0.8,
            by = 0.2
        ),
        labels = sprintf(
            "%.1f",
            seq(
                0,
                0.8,
                by = 0.2
            )
        )
    )

    # Very light grid.
    abline(
        v = seq(
            0,
            0.8,
            by = 0.2
        ),
        h = seq(
            0,
            0.8,
            by = 0.2
        ),
        col = "grey92",
        lty = 3,
        lwd = 0.8
    )

    # y = x reference line.
    abline(
        a = 0,
        b = 1,
        lty = 2,
        lwd = 1.3,
        col = "grey45"
    )

    # Region labels.
    text(
        x = 0.57,
        y = 0.80,
        labels =
            expression(
                epsilon^2[celltype] >
                epsilon^2[replicate]
            ),
        cex = 0.82,
        col = "grey35"
    )

    text(
        x = 0.70,
        y = 0.12,
        labels =
            expression(
                epsilon^2[replicate] >
                epsilon^2[celltype]
            ),
        cex = 0.82,
        col = "grey55"
    )

    # All components except CELL03.
    idx_other <-
        comp$CELL_component !=
        "CELL03"

    points(
        comp$replicate_epsilon2[
            idx_other
        ],
        comp$celltype_epsilon2[
            idx_other
        ],
        pch = 21,
        bg = "grey85",
        col = "grey20",
        lwd = 1.1,
        cex = point_cex[
            idx_other
        ]
    )

    # Highlight CELL03.
    idx3 <-
        comp$CELL_component ==
        "CELL03"

    points(
        comp$replicate_epsilon2[
            idx3
        ],
        comp$celltype_epsilon2[
            idx3
        ],
        pch = 21,
        bg = "firebrick2",
        col = "black",
        lwd = 1.3,
        cex = point_cex[
            idx3
        ] * 1.05
    )

    # Labels.
    for (
        i in seq_len(
            nrow(
                comp
            )
        )
    ) {

        lab <- paste0(
            comp$CELL_component[i],
            " (",
            comp$n_top100_cores[i],
            ")"
        )

        text(
            x =
                comp$label_x[i],
            y =
                comp$label_y[i],
            labels =
                lab,
            adj =
                c(
                    0,
                    0.5
                ),
            cex =
                if (
                    comp$CELL_component[i] ==
                    "CELL03"
                ) {
                    0.84
                } else {
                    0.76
                },
            font =
                if (
                    comp$CELL_component[i] ==
                    "CELL03"
                ) {
                    2
                } else {
                    1
                },
            col =
                if (
                    comp$CELL_component[i] ==
                    "CELL03"
                ) {
                    "firebrick4"
                } else {
                    "black"
                }
        )
    }

    # Compact summary in unused lower-right region.
    summary_lines <- c(
        paste0(
            n_above,
            "/",
            n_components,
            " components: cell-type > replicate"
        ),
        "100/100 top-core entries",
        paste0(
            "one-sided sign test: P = ",
            format(
                sign_test$p.value,
                scientific = TRUE,
                digits = 3
            )
        )
    )

    text(
        x = 0.43,
        y = 0.055,
        labels =
            paste(
                summary_lines,
                collapse = "\n"
            ),
        adj =
            c(
                0,
                0
            ),
        cex = 0.78,
        font = 2
    )

    # Small note explaining numbers in parentheses.
    mtext(
        "Numbers in parentheses indicate frequency among the 100 largest core elements.",
        side = 1,
        line = 4.15,
        cex = 0.72,
        col = "grey30"
    )
}


# ------------------------------------------------------------
# 6. Save PDF
# ------------------------------------------------------------

pdf(
    OUT_PDF,
    width = 7.3,
    height = 6.6,
    useDingbats = FALSE
)

draw_publication_figure()

dev.off()


# ------------------------------------------------------------
# 7. Save high-resolution PNG
# ------------------------------------------------------------

png(
    OUT_PNG,
    width = 2100,
    height = 1900,
    res = 300
)

draw_publication_figure()

dev.off()


# ------------------------------------------------------------
# 8. Console
# ------------------------------------------------------------

cat(
    "\nPublication figure completed.\n"
)

cat(
    "Distinct CELL components : ",
    n_components,
    "\n",
    sep = ""
)

cat(
    "Cell-type > replicate    : ",
    n_above,
    "/",
    n_components,
    "\n",
    sep = ""
)

cat(
    "Sign-test P              : ",
    sign_test$p.value,
    "\n",
    sep = ""
)

cat(
    "\nSaved:\n",
    OUT_PDF,
    "\n",
    OUT_PNG,
    "\n",
    sep = ""
)
