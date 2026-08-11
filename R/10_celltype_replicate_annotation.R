# ============================================================
# 11_celltype_replicate_annotation.R
# GSE303006 / GENCODE M23
#
# Biological annotation of Stage-2 CELL components appearing
# in the top 100 non-K-mixed Tucker core entries.
#
# Main goals:
#   1) inspect GSE303006_charm_metadata_qcpass.tsv.gz columns
#   2) match metadata to the final 4,258 Stage-2 cells
#   3) identify a cell-type annotation column
#   4) test each selected CELL component against cell type
#   5) test each selected CELL component against replicate R1/R2/R3
#   6) find cell types enriched among the top 50 positive-loading
#      and top 50 negative-loading cells
#
# IMPORTANT:
#   The continuous association tests use ALL 4,258 cells.
#   Top-50 analyses are used only for interpretation/enrichment.
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 0. Files / settings
# ------------------------------------------------------------

META_FILE <-
    "GSE303006_charm_metadata_qcpass.tsv.gz"

TOP_CORE_FILE <-
    "GSE303006_M23_stage2_top100_core.tsv"

SMALL_FACTORS_FILE <-
    "GSE303006_M23_stage2_Tucker_small_factors.rds"

MANIFEST_FILE <-
    "GSE303006_M23_stage2_implicit_manifest.rds"

TOP_N_EACH_SIGN <- 50L

# If automatic detection chooses the wrong column, set manually,
# e.g. CELLTYPE_COL <- "celltype"
CELLTYPE_COL <- NULL

OUT_COLUMN_INFO <-
    "GSE303006_M23_metadata_column_inspection.tsv"

OUT_MATCHED_META <-
    "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_CELLTYPE_TEST <-
    "GSE303006_M23_CELLcomponent_celltype_KW.tsv"

OUT_CELLTYPE_SUMMARY <-
    "GSE303006_M23_CELLcomponent_celltype_summary.tsv.gz"

OUT_CELLTYPE_ENRICH <-
    "GSE303006_M23_CELLcomponent_celltype_top50_enrichment.tsv.gz"

OUT_REPLICATE_TEST <-
    "GSE303006_M23_CELLcomponent_replicate_KW.tsv"

OUT_REPLICATE_SUMMARY <-
    "GSE303006_M23_CELLcomponent_replicate_summary.tsv"

OUT_REPLICATE_ENRICH <-
    "GSE303006_M23_CELLcomponent_replicate_top50_enrichment.tsv"

OUT_CORE_ANNOTATED <-
    "GSE303006_M23_stage2_top100_core_cell_annotation.tsv"

OUT_RDS <-
    "GSE303006_M23_CELLcomponent_biological_annotation.rds"


# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

need <- c(
    META_FILE,
    TOP_CORE_FILE,
    SMALL_FACTORS_FILE,
    MANIFEST_FILE
)

missing <- need[
    !file.exists(need)
]

if (length(missing) > 0L) {
    stop(
        "Missing input file(s):\n",
        paste(missing, collapse = "\n")
    )
}


# ------------------------------------------------------------
# 2. Read metadata and inspect columns
# ------------------------------------------------------------

meta <- fread(
    META_FILE
)

cat(
    "Metadata dimensions :",
    nrow(meta),
    "x",
    ncol(meta),
    "\n"
)

cat(
    "\nMetadata columns:\n"
)

print(
    names(meta)
)

# Basic information for every metadata column.
column_info <- rbindlist(
    lapply(
        names(meta),
        function(nm) {

            x <- meta[[nm]]

            data.table(
                column = nm,
                class =
                    paste(
                        class(x),
                        collapse = ","
                    ),
                n_nonNA =
                    sum(
                        !is.na(x)
                    ),
                n_unique =
                    uniqueN(
                        x[
                            !is.na(x)
                        ]
                    ),
                example =
                    paste(
                        head(
                            unique(
                                as.character(
                                    x[
                                        !is.na(x)
                                    ]
                                )
                            ),
                            5
                        ),
                        collapse = " | "
                    )
            )
        }
    )
)

fwrite(
    column_info,
    OUT_COLUMN_INFO,
    sep = "\t"
)

cat(
    "\nColumn inspection:\n"
)

print(
    column_info
)


# ------------------------------------------------------------
# 3. Load Stage-2 cell basis / final cells / top cores
# ------------------------------------------------------------

F <- readRDS(
    SMALL_FACTORS_FILE
)

M <- readRDS(
    MANIFEST_FILE
)

top_core <- fread(
    TOP_CORE_FILE
)

U_CELL <- F$U_CELL

if (is.null(U_CELL)) {
    stop(
        "U_CELL is absent from small_factors."
    )
}

final_cells <- rownames(
    U_CELL
)

if (is.null(final_cells)) {
    stop(
        "U_CELL has no cell row names."
    )
}

if (!identical(
    final_cells,
    M$cells
)) {
    stop(
        "Cell order differs between U_CELL and manifest."
    )
}

cat(
    "\nFinal Stage-2 cells :",
    length(final_cells),
    "\n"
)

selected_b <- sort(
    unique(
        as.integer(
            top_core$b_CELL
        )
    )
)

cat(
    "CELL components represented in top 100 CORE entries :",
    length(selected_b),
    "\n"
)

cat(
    "Components :",
    paste(
        selected_b,
        collapse = ", "
    ),
    "\n"
)


# ------------------------------------------------------------
# 4. Detect metadata cell-ID column by actual overlap
#
# We do not guess from the column name alone.
# For each metadata column, count overlap with final cell IDs.
# ------------------------------------------------------------

cell_id_candidates <- rbindlist(
    lapply(
        names(meta),
        function(nm) {

            x <- as.character(
                meta[[nm]]
            )

            ov <- sum(
                unique(
                    x[
                        !is.na(x)
                    ]
                ) %in%
                    final_cells
            )

            data.table(
                column = nm,
                overlap_final_cells = ov
            )
        }
    )
)

setorder(
    cell_id_candidates,
    -overlap_final_cells
)

cat(
    "\nBest metadata cell-ID candidates:\n"
)

print(
    head(
        cell_id_candidates,
        10
    )
)

CELL_ID_COL <-
    cell_id_candidates$column[1]

best_overlap <-
    cell_id_candidates$overlap_final_cells[1]

if (best_overlap == 0L) {
    stop(
        "No metadata column overlaps the final Stage-2 cell IDs."
    )
}

cat(
    "\nSelected CELL_ID_COL :",
    CELL_ID_COL,
    "\n"
)

cat(
    "Matched final cell IDs :",
    best_overlap,
    "/",
    length(final_cells),
    "\n"
)


# ------------------------------------------------------------
# 5. Restrict metadata to final 4,258 cells
# ------------------------------------------------------------

meta[
    ,
    .cell_id_for_join :=
        as.character(
            get(
                CELL_ID_COL
            )
        )
]

if (anyDuplicated(
    meta$.cell_id_for_join
)) {

    dup <- meta[
        duplicated(
            .cell_id_for_join
        ) |
        duplicated(
            .cell_id_for_join,
            fromLast = TRUE
        )
    ]

    stop(
        "Metadata cell-ID column contains duplicated IDs. ",
        "Examples: ",
        paste(
            head(
                unique(
                    dup$.cell_id_for_join
                ),
                10
            ),
            collapse = ", "
        )
    )
}

meta_final <- meta[
    match(
        final_cells,
        .cell_id_for_join
    )
]

meta_final[
    ,
    cell :=
        final_cells
]

matched <- !is.na(
    meta_final$.cell_id_for_join
)

cat(
    "\nMetadata matched to final cells :",
    sum(matched),
    "/",
    length(final_cells),
    "\n"
)

if (sum(matched) < 0.9 * length(final_cells)) {
    warning(
        "Less than 90% of final cells matched metadata."
    )
}


# ------------------------------------------------------------
# 6. Derive replicate from cell IDs
#
# GEO defines R1 / R2 / R3 prefixes as brain replicates.
# ------------------------------------------------------------

meta_final[
    ,
    replicate :=
        fifelse(
            grepl(
                "^R1",
                cell
            ),
            "R1",
            fifelse(
                grepl(
                    "^R2",
                    cell
                ),
                "R2",
                fifelse(
                    grepl(
                        "^R3",
                        cell
                    ),
                    "R3",
                    NA_character_
                )
            )
        )
]

cat(
    "\nReplicate counts in final cells:\n"
)

print(
    table(
        meta_final$replicate,
        useNA = "ifany"
    )
)


# ------------------------------------------------------------
# 7. Detect likely cell-type column
#
# Strategy:
#   - prefer names containing celltype / cell_type / annotation /
#     subtype / cluster / class
#   - require a categorical-looking number of values
#   - score by column name + coverage
#
# If this chooses incorrectly, set CELLTYPE_COL manually above.
# ------------------------------------------------------------

if (is.null(
    CELLTYPE_COL
)) {

    patt <-
        "cell.?type|celltype|annotation|annot|subtype|cluster|cell.?class|class|identity|ident"

    candidate_names <- names(meta_final)[
        grepl(
            patt,
            names(meta_final),
            ignore.case = TRUE
        )
    ]

    candidate_names <- setdiff(
        candidate_names,
        c(
            CELL_ID_COL,
            ".cell_id_for_join",
            "cell",
            "replicate"
        )
    )

    if (length(candidate_names) == 0L) {

        # Fallback: categorical columns with 2-100 unique levels.
        candidate_names <- names(meta_final)[
            vapply(
                meta_final,
                function(x) {

                    nlev <- uniqueN(
                        x[
                            !is.na(x)
                        ]
                    )

                    nlev >= 2L &&
                        nlev <= 100L
                },
                logical(1)
            )
        ]

        candidate_names <- setdiff(
            candidate_names,
            c(
                CELL_ID_COL,
                ".cell_id_for_join",
                "cell",
                "replicate"
            )
        )
    }

    if (length(candidate_names) == 0L) {
        stop(
            "No plausible cell-type annotation column was detected. ",
            "Inspect ",
            OUT_COLUMN_INFO,
            " and set CELLTYPE_COL manually."
        )
    }

    celltype_candidates <- rbindlist(
        lapply(
            candidate_names,
            function(nm) {

                x <- meta_final[[nm]]

                n_nonNA <- sum(
                    !is.na(x)
                )

                n_unique <- uniqueN(
                    x[
                        !is.na(x)
                    ]
                )

                name_score <-
                    10 *
                    grepl(
                        "cell.?type|celltype",
                        nm,
                        ignore.case = TRUE
                    ) +
                    5 *
                    grepl(
                        "annotation|annot|subtype",
                        nm,
                        ignore.case = TRUE
                    ) +
                    2 *
                    grepl(
                        "cluster|class|identity|ident",
                        nm,
                        ignore.case = TRUE
                    )

                coverage_score <-
                    n_nonNA /
                    nrow(meta_final)

                cardinality_score <-
                    as.numeric(
                        n_unique >= 2L &&
                        n_unique <= 100L
                    )

                data.table(
                    column = nm,
                    n_nonNA = n_nonNA,
                    n_unique = n_unique,
                    score =
                        name_score +
                        coverage_score +
                        cardinality_score
                )
            }
        )
    )

    setorder(
        celltype_candidates,
        -score,
        n_unique
    )

    cat(
        "\nCell-type annotation candidates:\n"
    )

    print(
        celltype_candidates
    )

    CELLTYPE_COL <-
        celltype_candidates$column[1]
}

if (!CELLTYPE_COL %in%
    names(meta_final)) {
    stop(
        "CELLTYPE_COL not found: ",
        CELLTYPE_COL
    )
}

meta_final[
    ,
    cell_type :=
        as.character(
            get(
                CELLTYPE_COL
            )
        )
]

cat(
    "\nSelected CELLTYPE_COL :",
    CELLTYPE_COL,
    "\n"
)

cat(
    "Number of annotated cell types :",
    uniqueN(
        meta_final$cell_type[
            !is.na(
                meta_final$cell_type
            )
        ]
    ),
    "\n"
)

cat(
    "\nCell-type counts:\n"
)

print(
    sort(
        table(
            meta_final$cell_type,
            useNA = "ifany"
        ),
        decreasing = TRUE
    )
)

fwrite(
    meta_final,
    OUT_MATCHED_META,
    sep = "\t"
)


# ------------------------------------------------------------
# 8. Helper functions
# ------------------------------------------------------------

kw_epsilon2 <- function(
    x,
    group
) {

    ok <-
        is.finite(x) &
        !is.na(group)

    x <- x[ok]
    group <- factor(
        group[ok]
    )

    if (length(x) < 2L ||
        nlevels(group) < 2L) {
        return(
            list(
                statistic = NA_real_,
                df = NA_real_,
                p = NA_real_,
                epsilon2 = NA_real_
            )
        )
    }

    kt <- kruskal.test(
        x ~ group
    )

    H <- as.numeric(
        kt$statistic
    )

    klev <- nlevels(
        group
    )

    n <- length(
        x
    )

    eps2 <- (
        H -
        klev +
        1
    ) /
        (
            n -
            klev
        )

    eps2 <- max(
        0,
        min(
            1,
            eps2
        )
    )

    list(
        statistic = H,
        df = as.numeric(
            kt$parameter
        ),
        p = kt$p.value,
        epsilon2 = eps2
    )
}


fisher_one_level <- function(
    selected,
    group,
    level
) {

    ok <- !is.na(group)

    selected <- selected[ok]
    group <- group[ok]

    a <- sum(
        selected &
        group == level
    )

    b <- sum(
        selected &
        group != level
    )

    c <- sum(
        !selected &
        group == level
    )

    d <- sum(
        !selected &
        group != level
    )

    mat <- matrix(
        c(
            a,
            b,
            c,
            d
        ),
        nrow = 2,
        byrow = TRUE
    )

    ft <- fisher.test(
        mat,
        alternative = "greater"
    )

    data.table(
        selected_in_level = a,
        selected_not_level = b,
        background_in_level = c,
        background_not_level = d,
        odds_ratio =
            unname(
                ft$estimate
            ),
        p_value =
            ft$p.value
    )
}


# ------------------------------------------------------------
# 9. Continuous cell-type association for selected CELL modes
# ------------------------------------------------------------

celltype_tests <- list()
celltype_summaries <- list()

replicate_tests <- list()
replicate_summaries <- list()

ctt <- 1L
cts <- 1L
rtt <- 1L
rts <- 1L

for (b in selected_b) {

    loading <- U_CELL[
        ,
        b
    ]

    # Cell type: all cells
    kw_ct <- kw_epsilon2(
        loading,
        meta_final$cell_type
    )

    celltype_tests[[ctt]] <- data.table(
        b_CELL = b,
        CELL_component =
            sprintf(
                "CELL%02d",
                b
            ),
        n_cells =
            sum(
                is.finite(loading) &
                !is.na(
                    meta_final$cell_type
                )
            ),
        n_celltypes =
            uniqueN(
                meta_final$cell_type[
                    !is.na(
                        meta_final$cell_type
                    )
                ]
            ),
        KW_statistic =
            kw_ct$statistic,
        KW_df =
            kw_ct$df,
        p_value =
            kw_ct$p,
        epsilon2 =
            kw_ct$epsilon2
    )

    ctt <- ctt + 1L

    tmp_ct <- data.table(
        cell_type =
            meta_final$cell_type,
        loading =
            loading
    )[
        !is.na(
            cell_type
        )
    ][
        ,
        .(
            n = .N,
            mean_loading =
                mean(
                    loading
                ),
            median_loading =
                median(
                    loading
                ),
            sd_loading =
                sd(
                    loading
                ),
            mean_abs_loading =
                mean(
                    abs(
                        loading
                    )
                )
        ),
        by = cell_type
    ]

    tmp_ct[
        ,
        b_CELL := b
    ]

    tmp_ct[
        ,
        CELL_component :=
            sprintf(
                "CELL%02d",
                b
            )
    ]

    setcolorder(
        tmp_ct,
        c(
            "b_CELL",
            "CELL_component",
            "cell_type",
            "n",
            "mean_loading",
            "median_loading",
            "sd_loading",
            "mean_abs_loading"
        )
    )

    celltype_summaries[[cts]] <-
        tmp_ct

    cts <- cts + 1L

    # Replicate: all cells
    kw_rep <- kw_epsilon2(
        loading,
        meta_final$replicate
    )

    replicate_tests[[rtt]] <- data.table(
        b_CELL = b,
        CELL_component =
            sprintf(
                "CELL%02d",
                b
            ),
        n_cells =
            sum(
                is.finite(loading) &
                !is.na(
                    meta_final$replicate
                )
            ),
        n_replicates =
            uniqueN(
                meta_final$replicate[
                    !is.na(
                        meta_final$replicate
                    )
                ]
            ),
        KW_statistic =
            kw_rep$statistic,
        KW_df =
            kw_rep$df,
        p_value =
            kw_rep$p,
        epsilon2 =
            kw_rep$epsilon2
    )

    rtt <- rtt + 1L

    tmp_rep <- data.table(
        replicate =
            meta_final$replicate,
        loading =
            loading
    )[
        !is.na(
            replicate
        )
    ][
        ,
        .(
            n = .N,
            mean_loading =
                mean(
                    loading
                ),
            median_loading =
                median(
                    loading
                ),
            sd_loading =
                sd(
                    loading
                ),
            mean_abs_loading =
                mean(
                    abs(
                        loading
                    )
                )
        ),
        by = replicate
    ]

    tmp_rep[
        ,
        b_CELL := b
    ]

    tmp_rep[
        ,
        CELL_component :=
            sprintf(
                "CELL%02d",
                b
            )
    ]

    setcolorder(
        tmp_rep,
        c(
            "b_CELL",
            "CELL_component",
            "replicate",
            "n",
            "mean_loading",
            "median_loading",
            "sd_loading",
            "mean_abs_loading"
        )
    )

    replicate_summaries[[rts]] <-
        tmp_rep

    rts <- rts + 1L
}

celltype_test <- rbindlist(
    celltype_tests
)

celltype_test[
    ,
    FDR :=
        p.adjust(
            p_value,
            method = "BH"
        )
]

setorder(
    celltype_test,
    p_value
)

celltype_summary <- rbindlist(
    celltype_summaries
)

replicate_test <- rbindlist(
    replicate_tests
)

replicate_test[
    ,
    FDR :=
        p.adjust(
            p_value,
            method = "BH"
        )
]

setorder(
    replicate_test,
    p_value
)

replicate_summary <- rbindlist(
    replicate_summaries
)

fwrite(
    celltype_test,
    OUT_CELLTYPE_TEST,
    sep = "\t"
)

fwrite(
    celltype_summary,
    OUT_CELLTYPE_SUMMARY,
    sep = "\t"
)

fwrite(
    replicate_test,
    OUT_REPLICATE_TEST,
    sep = "\t"
)

fwrite(
    replicate_summary,
    OUT_REPLICATE_SUMMARY,
    sep = "\t"
)


# ------------------------------------------------------------
# 10. Top-50 positive / negative cell-type enrichment
#
# For each CELL component, select:
#   top 50 positive-loading cells
#   top 50 negative-loading cells
#
# Fisher exact test compares each cell type with all other cells.
# ------------------------------------------------------------

ct_enrich_list <- list()
rep_enrich_list <- list()

cte <- 1L
rpe <- 1L

for (b in selected_b) {

    loading <- U_CELL[
        ,
        b
    ]

    n_take <- min(
        TOP_N_EACH_SIGN,
        length(loading)
    )

    pos_idx <- order(
        loading,
        decreasing = TRUE
    )[
        seq_len(
            n_take
        )
    ]

    neg_idx <- order(
        loading,
        decreasing = FALSE
    )[
        seq_len(
            n_take
        )
    ]

    for (direction in c(
        "positive",
        "negative"
    )) {

        idx <- if (
            direction == "positive"
        ) {
            pos_idx
        } else {
            neg_idx
        }

        selected <- rep(
            FALSE,
            length(loading)
        )

        selected[idx] <- TRUE

        # Cell type enrichment
        levels_ct <- sort(
            unique(
                meta_final$cell_type[
                    !is.na(
                        meta_final$cell_type
                    )
                ]
            )
        )

        for (lev in levels_ct) {

            ft <- fisher_one_level(
                selected,
                meta_final$cell_type,
                lev
            )

            ft[
                ,
                b_CELL := b
            ]

            ft[
                ,
                CELL_component :=
                    sprintf(
                        "CELL%02d",
                        b
                    )
            ]

            ft[
                ,
                direction :=
                    direction
            ]

            ft[
                ,
                cell_type :=
                    lev
            ]

            selected_total <-
                sum(
                    selected &
                    !is.na(
                        meta_final$cell_type
                    )
                )

            ft[
                ,
                selected_total :=
                    selected_total
            ]

            ft[
                ,
                fraction_selected :=
                    if (
                        selected_total > 0
                    ) {
                        selected_in_level /
                            selected_total
                    } else {
                        NA_real_
                    }
            ]

            ct_enrich_list[[cte]] <-
                ft

            cte <- cte + 1L
        }

        # Replicate enrichment
        levels_rep <- c(
            "R1",
            "R2",
            "R3"
        )

        for (lev in levels_rep) {

            ft <- fisher_one_level(
                selected,
                meta_final$replicate,
                lev
            )

            ft[
                ,
                b_CELL := b
            ]

            ft[
                ,
                CELL_component :=
                    sprintf(
                        "CELL%02d",
                        b
                    )
            ]

            ft[
                ,
                direction :=
                    direction
            ]

            ft[
                ,
                replicate :=
                    lev
            ]

            selected_total <-
                sum(
                    selected &
                    !is.na(
                        meta_final$replicate
                    )
                )

            ft[
                ,
                selected_total :=
                    selected_total
            ]

            ft[
                ,
                fraction_selected :=
                    if (
                        selected_total > 0
                    ) {
                        selected_in_level /
                            selected_total
                    } else {
                        NA_real_
                    }
            ]

            rep_enrich_list[[rpe]] <-
                ft

            rpe <- rpe + 1L
        }
    }
}

celltype_enrich <- rbindlist(
    ct_enrich_list
)

celltype_enrich[
    ,
    FDR :=
        p.adjust(
            p_value,
            method = "BH"
        )
]

setorder(
    celltype_enrich,
    b_CELL,
    direction,
    p_value
)

replicate_enrich <- rbindlist(
    rep_enrich_list
)

replicate_enrich[
    ,
    FDR :=
        p.adjust(
            p_value,
            method = "BH"
        )
]

setorder(
    replicate_enrich,
    b_CELL,
    direction,
    p_value
)

fwrite(
    celltype_enrich,
    OUT_CELLTYPE_ENRICH,
    sep = "\t"
)

fwrite(
    replicate_enrich,
    OUT_REPLICATE_ENRICH,
    sep = "\t"
)


# ------------------------------------------------------------
# 11. Add a compact biological annotation to top 100 cores
#
# For each b_CELL:
#   - cell type with largest mean positive loading
#   - cell type with most negative mean loading
#   - strongest positive/negative top-50 enrichment
#   - cell-type KW effect size
#   - replicate KW effect size
# ------------------------------------------------------------

best_mean_pos <- celltype_summary[
    order(
        b_CELL,
        -mean_loading
    ),
    .SD[1],
    by = b_CELL
][
    ,
    .(
        b_CELL,
        celltype_high_mean =
            cell_type,
        celltype_high_mean_loading =
            mean_loading
    )
]

best_mean_neg <- celltype_summary[
    order(
        b_CELL,
        mean_loading
    ),
    .SD[1],
    by = b_CELL
][
    ,
    .(
        b_CELL,
        celltype_low_mean =
            cell_type,
        celltype_low_mean_loading =
            mean_loading
    )
]

best_enrich_pos <- celltype_enrich[
    direction == "positive"
][
    order(
        b_CELL,
        p_value,
        -odds_ratio
    ),
    .SD[1],
    by = b_CELL
][
    ,
    .(
        b_CELL,
        positive_top50_celltype =
            cell_type,
        positive_top50_OR =
            odds_ratio,
        positive_top50_p =
            p_value,
        positive_top50_FDR =
            FDR,
        positive_top50_fraction =
            fraction_selected
    )
]

best_enrich_neg <- celltype_enrich[
    direction == "negative"
][
    order(
        b_CELL,
        p_value,
        -odds_ratio
    ),
    .SD[1],
    by = b_CELL
][
    ,
    .(
        b_CELL,
        negative_top50_celltype =
            cell_type,
        negative_top50_OR =
            odds_ratio,
        negative_top50_p =
            p_value,
        negative_top50_FDR =
            FDR,
        negative_top50_fraction =
            fraction_selected
    )
]

ct_test_compact <- celltype_test[
    ,
    .(
        b_CELL,
        celltype_KW_p =
            p_value,
        celltype_KW_FDR =
            FDR,
        celltype_KW_epsilon2 =
            epsilon2
    )
]

rep_test_compact <- replicate_test[
    ,
    .(
        b_CELL,
        replicate_KW_p =
            p_value,
        replicate_KW_FDR =
            FDR,
        replicate_KW_epsilon2 =
            epsilon2
    )
]

core_annotated <- merge(
    copy(top_core),
    best_mean_pos,
    by = "b_CELL",
    all.x = TRUE,
    sort = FALSE
)

core_annotated <- merge(
    core_annotated,
    best_mean_neg,
    by = "b_CELL",
    all.x = TRUE,
    sort = FALSE
)

core_annotated <- merge(
    core_annotated,
    best_enrich_pos,
    by = "b_CELL",
    all.x = TRUE,
    sort = FALSE
)

core_annotated <- merge(
    core_annotated,
    best_enrich_neg,
    by = "b_CELL",
    all.x = TRUE,
    sort = FALSE
)

core_annotated <- merge(
    core_annotated,
    ct_test_compact,
    by = "b_CELL",
    all.x = TRUE,
    sort = FALSE
)

core_annotated <- merge(
    core_annotated,
    rep_test_compact,
    by = "b_CELL",
    all.x = TRUE,
    sort = FALSE
)

setorder(
    core_annotated,
    core_rank
)

fwrite(
    core_annotated,
    OUT_CORE_ANNOTATED,
    sep = "\t"
)


# ------------------------------------------------------------
# 12. Save RDS
# ------------------------------------------------------------

result <- list(
    metadata_columns =
        column_info,
    cell_id_column =
        CELL_ID_COL,
    celltype_column =
        CELLTYPE_COL,
    metadata_final =
        meta_final,
    selected_CELL_components =
        selected_b,
    celltype_test =
        celltype_test,
    celltype_summary =
        celltype_summary,
    celltype_top50_enrichment =
        celltype_enrich,
    replicate_test =
        replicate_test,
    replicate_summary =
        replicate_summary,
    replicate_top50_enrichment =
        replicate_enrich,
    top100_core_annotated =
        core_annotated
)

saveRDS(
    result,
    OUT_RDS
)


# ------------------------------------------------------------
# 13. Console summary
# ------------------------------------------------------------

cat(
    "\n============================================================\n"
)

cat(
    "CELL component biological annotation completed\n"
)

cat(
    "============================================================\n"
)

cat(
    "Metadata cell ID column :",
    CELL_ID_COL,
    "\n"
)

cat(
    "Metadata cell-type column :",
    CELLTYPE_COL,
    "\n"
)

cat(
    "Matched final cells :",
    sum(matched),
    "/",
    length(final_cells),
    "\n"
)

cat(
    "\nCell-type association tests:\n"
)

print(
    celltype_test[
        order(
            p_value
        )
    ]
)

cat(
    "\nReplicate association tests:\n"
)

print(
    replicate_test[
        order(
            p_value
        )
    ]
)

cat(
    "\nMost interpretable top-20 CORE entries:\n"
)

print(
    core_annotated[
        core_rank <= 20
    ][
        ,
        .(
            core_rank,
            b_CELL,
            core_value,
            positive_top50_celltype,
            positive_top50_OR,
            positive_top50_FDR,
            negative_top50_celltype,
            negative_top50_OR,
            negative_top50_FDR,
            celltype_KW_epsilon2,
            celltype_KW_FDR,
            replicate_KW_epsilon2,
            replicate_KW_FDR
        )
    ]
)

cat(
    "\nSaved files:\n",
    OUT_COLUMN_INFO, "\n",
    OUT_MATCHED_META, "\n",
    OUT_CELLTYPE_TEST, "\n",
    OUT_CELLTYPE_SUMMARY, "\n",
    OUT_CELLTYPE_ENRICH, "\n",
    OUT_REPLICATE_TEST, "\n",
    OUT_REPLICATE_SUMMARY, "\n",
    OUT_REPLICATE_ENRICH, "\n",
    OUT_CORE_ANNOTATED, "\n",
    OUT_RDS, "\n",
    sep = ""
)
