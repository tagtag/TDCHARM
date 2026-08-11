# ============================================================
# 15_core12_GO_BP_reduce_redundancy_leading_edge.R
# GSE303006 / GENCODE M23
#
# Reduce redundancy among significant Mouse GO Biological
# Process GSEA terms for core #12 and report representative
# biological processes with their leading-edge genes.
#
# Strategy:
#   - Positive and negative NES terms are handled separately.
#   - Only FDR <= 0.05 terms are considered.
#   - GO semantic similarity is calculated with GOSemSim,
#     Wang method, BP ontology.
#   - Terms with semantic similarity >= 0.70 are grouped.
#   - Greedy grouping is performed after ordering terms by:
#         1) smaller FDR
#         2) larger |NES|
#     Therefore the first term in each group is the
#     representative term.
#   - Up to 10 representative processes are reported per side.
#   - Leading-edge genes are taken directly from the GSEA
#     core_enrichment field.
#
# Outputs include:
#   1) all significant terms with redundancy-group assignment
#   2) all nonredundant representative terms
#   3) top 10 representative terms for positive and negative sides
#   4) long-format leading-edge genes for those representatives
# ============================================================


# ============================================================
# 0. Settings
# ============================================================

GSEA_RDS <-
    "GSE303006_M23_core12_GO_BP_GSEA.rds"

RNK_FILE <-
    "GSE303006_M23_core12_all16239_gene_rank_degree_corrected.rnk"

FDR_CUTOFF <- 0.05

SEMANTIC_CUTOFF <- 0.70

TOP_REPRESENTATIVE <- 10L

SEMANTIC_METHOD <- "Wang"


# Output files
OUT_POS_GROUPED <-
    "GSE303006_M23_core12_GO_BP_positive_FDR005_redundancy_groups.tsv"

OUT_NEG_GROUPED <-
    "GSE303006_M23_core12_GO_BP_negative_FDR005_redundancy_groups.tsv"

OUT_POS_REP_ALL <-
    "GSE303006_M23_core12_GO_BP_positive_nonredundant_all.tsv"

OUT_NEG_REP_ALL <-
    "GSE303006_M23_core12_GO_BP_negative_nonredundant_all.tsv"

OUT_POS_TOP <-
    "GSE303006_M23_core12_GO_BP_positive_representative_top10.tsv"

OUT_NEG_TOP <-
    "GSE303006_M23_core12_GO_BP_negative_representative_top10.tsv"

OUT_POS_LEADING_LONG <-
    "GSE303006_M23_core12_GO_BP_positive_top10_leading_edge_genes.tsv"

OUT_NEG_LEADING_LONG <-
    "GSE303006_M23_core12_GO_BP_negative_top10_leading_edge_genes.tsv"

OUT_POS_SIM <-
    "GSE303006_M23_core12_GO_BP_positive_semantic_similarity.rds"

OUT_NEG_SIM <-
    "GSE303006_M23_core12_GO_BP_negative_semantic_similarity.rds"

OUT_RDS <-
    "GSE303006_M23_core12_GO_BP_redundancy_reduced.rds"


# ============================================================
# 1. Package checks
# ============================================================

needed <- c(
    "data.table",
    "GOSemSim",
    "org.Mm.eg.db",
    "AnnotationDbi",
    "GO.db"
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
        paste0(
            "Missing package(s): ",
            paste(missing, collapse = ", "),
            "\n\nInstall with:\n",
            "if (!requireNamespace(\"BiocManager\", quietly=TRUE)) ",
            "install.packages(\"BiocManager\")\n",
            "BiocManager::install(c(",
            paste(
                sprintf("\"%s\"", setdiff(missing, "data.table")),
                collapse = ", "
            ),
            "))\n",
            if ("data.table" %in% missing) {
                "install.packages(\"data.table\")\n"
            } else {
                ""
            }
        )
    )
}

library(data.table)


# ============================================================
# 2. Input checks
# ============================================================

need_files <- c(
    GSEA_RDS,
    RNK_FILE
)

missing_files <- need_files[
    !file.exists(need_files)
]

if (length(missing_files) > 0L) {
    stop(
        "Missing input file(s):\n",
        paste(missing_files, collapse = "\n")
    )
}


# ============================================================
# 3. Load GSEA results
# ============================================================

gsea_bundle <- readRDS(
    GSEA_RDS
)

if (!"all_results" %in% names(gsea_bundle)) {
    stop(
        "GSEA RDS does not contain all_results."
    )
}

res <- as.data.table(
    copy(
        gsea_bundle$all_results
    )
)

required_cols <- c(
    "ID",
    "Description",
    "setSize",
    "NES",
    "pvalue",
    "p.adjust",
    "core_enrichment"
)

if (!all(required_cols %in% names(res))) {
    stop(
        "Missing GSEA result columns: ",
        paste(
            setdiff(required_cols, names(res)),
            collapse = ", "
        )
    )
}

pos <- res[
    NES > 0 &
    p.adjust <= FDR_CUTOFF
]

neg <- res[
    NES < 0 &
    p.adjust <= FDR_CUTOFF
]

cat(
    "Significant positive GO BP terms :",
    nrow(pos),
    "\n"
)

cat(
    "Significant negative GO BP terms :",
    nrow(neg),
    "\n"
)

if (nrow(pos) == 0L || nrow(neg) == 0L) {
    stop(
        "One side contains no significant GO BP terms."
    )
}


# ============================================================
# 4. Load ranked genes
# ============================================================

rank_dt <- fread(
    RNK_FILE,
    header = FALSE,
    col.names = c(
        "gene_name",
        "rank_score"
    )
)

rank_dt[
    ,
    gene_name :=
        as.character(
            gene_name
        )
]

rank_dt[
    ,
    rank_score :=
        as.numeric(
            rank_score
        )
]

setorder(
    rank_dt,
    -rank_score,
    gene_name
)

rank_dt[
    ,
    global_signed_rank :=
        seq_len(.N)
]

score_lookup <- rank_dt$rank_score
names(score_lookup) <- rank_dt$gene_name

rank_lookup <- rank_dt$global_signed_rank
names(rank_lookup) <- rank_dt$gene_name


# ============================================================
# 5. Prepare semantic-data object
#
# Wang similarity is graph-based, so computeIC=FALSE is enough.
# ============================================================

cat(
    "\nPreparing Mouse GO BP semantic data...\n"
)

# Use the first positional argument for maximum compatibility
# across older and newer GOSemSim versions.
# Older GOSemSim versions call this argument OrgDb; newer versions
# also accept it positionally (and may issue only a deprecation warning).
semData <- GOSemSim::godata(
    "org.Mm.eg.db",
    ont =
        "BP",
    computeIC =
        FALSE
)


# ============================================================
# 6. Helper: calculate GO semantic similarity matrix
# ============================================================

make_semantic_matrix <- function(
    ids
) {

    ids <- unique(
        as.character(
            ids
        )
    )

    sim <- GOSemSim::mgoSim(
        GO1 =
            ids,
        GO2 =
            ids,
        semData =
            semData,
        measure =
            SEMANTIC_METHOD,
        combine =
            NULL
    )

    sim <- as.matrix(
        sim
    )

    # Force expected row/column order if possible.
    sim <- sim[
        ids,
        ids,
        drop = FALSE
    ]

    # Invalid/unknown similarities should not trigger grouping.
    sim[
        !is.finite(sim)
    ] <- 0

    diag(
        sim
    ) <- 1

    sim
}


# ============================================================
# 7. Helper: greedy semantic redundancy grouping
#
# Terms are ordered by:
#   smallest FDR, then largest |NES|.
#
# For each term:
#   - compare with all existing representative terms
#   - if maximum similarity >= cutoff, assign it to the
#     most similar representative group
#   - otherwise create a new group and make this term the
#     representative
#
# Thus every representative is the most statistically
# significant term encountered for its semantic neighborhood.
# ============================================================

reduce_one_side <- function(
    dt,
    direction
) {

    x <- copy(
        dt
    )

    x[
        ,
        abs_NES :=
            abs(
                NES
            )
    ]

    setorder(
        x,
        p.adjust,
        -abs_NES,
        ID
    )

    ids <- x$ID

    cat(
        "\nCalculating semantic similarity for ",
        direction,
        " terms (n=",
        length(ids),
        ")...\n",
        sep = ""
    )

    sim <- make_semantic_matrix(
        ids
    )

    representative_ids <- character(0)

    group_id <- integer(
        nrow(x)
    )

    representative_id_for_term <- character(
        nrow(x)
    )

    max_similarity_to_rep <- numeric(
        nrow(x)
    )

    for (i in seq_len(
        nrow(x)
    )) {

        this_id <- x$ID[i]

        if (length(
            representative_ids
        ) == 0L) {

            representative_ids <-
                this_id

            group_id[i] <- 1L

            representative_id_for_term[i] <-
                this_id

            max_similarity_to_rep[i] <-
                1

            next
        }

        ss <- sim[
            this_id,
            representative_ids,
            drop = TRUE
        ]

        best_j <- which.max(
            ss
        )

        best_sim <- ss[
            best_j
        ]

        if (
            length(best_sim) == 1L &&
            is.finite(best_sim) &&
            best_sim >= SEMANTIC_CUTOFF
        ) {

            group_id[i] <-
                best_j

            representative_id_for_term[i] <-
                representative_ids[
                    best_j
                ]

            max_similarity_to_rep[i] <-
                best_sim

        } else {

            representative_ids <-
                c(
                    representative_ids,
                    this_id
                )

            new_group <-
                length(
                    representative_ids
                )

            group_id[i] <-
                new_group

            representative_id_for_term[i] <-
                this_id

            max_similarity_to_rep[i] <-
                1
        }
    }

    x[
        ,
        redundancy_group :=
            group_id
    ]

    x[
        ,
        representative_GO :=
            representative_id_for_term
    ]

    x[
        ,
        semantic_similarity_to_representative :=
            max_similarity_to_rep
    ]

    x[
        ,
        is_representative :=
            ID ==
            representative_GO
    ]

    # Attach representative description.
    rep_desc <- x[
        is_representative == TRUE,
        .(
            representative_GO =
                ID,
            representative_Description =
                Description
        )
    ]

    x[
        rep_desc,
        on = "representative_GO",
        representative_Description :=
            i.representative_Description
    ]

    # Group size.
    x[
        ,
        redundancy_group_size :=
            .N,
        by = redundancy_group
    ]

    reps <- x[
        is_representative == TRUE
    ]

    # Number of terms represented by this representative.
    group_counts <- x[
        ,
        .(
            n_terms_in_semantic_group =
                .N
        ),
        by = redundancy_group
    ]

    reps[
        group_counts,
        on = "redundancy_group",
        n_terms_in_semantic_group :=
            i.n_terms_in_semantic_group
    ]

    setorder(
        reps,
        p.adjust,
        -abs_NES,
        ID
    )

    reps[
        ,
        representative_rank :=
            seq_len(.N)
    ]

    list(
        grouped =
            x,

        representatives =
            reps,

        similarity =
            sim
    )
}


# ============================================================
# 8. Perform redundancy reduction separately by sign
# ============================================================

pos_red <- reduce_one_side(
    pos,
    "positive"
)

neg_red <- reduce_one_side(
    neg,
    "negative"
)

saveRDS(
    pos_red$similarity,
    OUT_POS_SIM
)

saveRDS(
    neg_red$similarity,
    OUT_NEG_SIM
)

fwrite(
    pos_red$grouped,
    OUT_POS_GROUPED,
    sep = "\t"
)

fwrite(
    neg_red$grouped,
    OUT_NEG_GROUPED,
    sep = "\t"
)

fwrite(
    pos_red$representatives,
    OUT_POS_REP_ALL,
    sep = "\t"
)

fwrite(
    neg_red$representatives,
    OUT_NEG_REP_ALL,
    sep = "\t"
)


# ============================================================
# 9. Choose top representative processes
# ============================================================

pos_top <- pos_red$representatives[
    seq_len(
        min(
            TOP_REPRESENTATIVE,
            .N
        )
    )
]

neg_top <- neg_red$representatives[
    seq_len(
        min(
            TOP_REPRESENTATIVE,
            .N
        )
    )
]


# ============================================================
# 10. Helper: parse and rank leading-edge genes
#
# GSEA core_enrichment uses "/"-separated gene symbols.
# For positive terms, genes are ordered from highest positive
# core #12 ranking score downward.
# For negative terms, genes are ordered from the most negative
# ranking score upward.
# ============================================================

leading_edge_long <- function(
    top_dt,
    direction
) {

    ans <- vector(
        "list",
        nrow(top_dt)
    )

    for (i in seq_len(
        nrow(top_dt)
    )) {

        genes <- strsplit(
            as.character(
                top_dt$core_enrichment[i]
            ),
            "/",
            fixed = TRUE
        )[[1]]

        genes <- unique(
            genes[
                nzchar(
                    genes
                )
            ]
        )

        scores <- unname(
            score_lookup[
                genes
            ]
        )

        global_ranks <- unname(
            rank_lookup[
                genes
            ]
        )

        if (
            direction ==
            "positive"
        ) {

            oo <- order(
                scores,
                decreasing = TRUE,
                na.last = TRUE
            )

        } else {

            oo <- order(
                scores,
                decreasing = FALSE,
                na.last = TRUE
            )
        }

        genes <-
            genes[oo]

        scores <-
            scores[oo]

        global_ranks <-
            global_ranks[oo]

        ans[[i]] <- data.table(
            direction =
                direction,

            representative_rank =
                top_dt$representative_rank[i],

            GO_ID =
                top_dt$ID[i],

            Description =
                top_dt$Description[i],

            NES =
                top_dt$NES[i],

            p_adjust =
                top_dt$p.adjust[i],

            n_terms_in_semantic_group =
                top_dt$n_terms_in_semantic_group[i],

            leading_edge_rank =
                seq_along(
                    genes
                ),

            gene_name =
                genes,

            core12_gene_rank_score =
                scores,

            global_signed_rank =
                global_ranks
        )
    }

    rbindlist(
        ans,
        fill = TRUE
    )
}


pos_leading <- leading_edge_long(
    pos_top,
    "positive"
)

neg_leading <- leading_edge_long(
    neg_top,
    "negative"
)


# ============================================================
# 11. Add leading-edge summary columns to representative table
# ============================================================

make_rep_summary <- function(
    top_dt,
    leading_dt,
    direction
) {

    y <- copy(
        top_dt
    )

    lead_summary <- leading_dt[
        ,
        .(
            n_leading_edge_genes =
                .N,

            leading_edge_genes =
                paste(
                    gene_name,
                    collapse = "/"
                ),

            top10_leading_edge_genes =
                paste(
                    head(
                        gene_name,
                        10
                    ),
                    collapse = "/"
                )
        ),
        by = .(
            representative_rank,
            GO_ID
        )
    ]

    y[
        lead_summary,
        on = c(
            "representative_rank",
            "ID" =
                "GO_ID"
        ),
        `:=`(
            n_leading_edge_genes =
                i.n_leading_edge_genes,

            leading_edge_genes =
                i.leading_edge_genes,

            top10_leading_edge_genes =
                i.top10_leading_edge_genes
        )
    ]

    y[
        ,
        direction :=
            direction
    ]

    keep <- c(
        "direction",
        "representative_rank",
        "ID",
        "Description",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust",
        "n_terms_in_semantic_group",
        "n_leading_edge_genes",
        "top10_leading_edge_genes",
        "leading_edge_genes"
    )

    y[
        ,
        ..keep
    ]
}


pos_summary <- make_rep_summary(
    pos_top,
    pos_leading,
    "positive"
)

neg_summary <- make_rep_summary(
    neg_top,
    neg_leading,
    "negative"
)

fwrite(
    pos_summary,
    OUT_POS_TOP,
    sep = "\t"
)

fwrite(
    neg_summary,
    OUT_NEG_TOP,
    sep = "\t"
)

fwrite(
    pos_leading,
    OUT_POS_LEADING_LONG,
    sep = "\t"
)

fwrite(
    neg_leading,
    OUT_NEG_LEADING_LONG,
    sep = "\t"
)


# ============================================================
# 12. Console display
# ============================================================

show_cols <- c(
    "representative_rank",
    "ID",
    "Description",
    "NES",
    "p.adjust",
    "n_terms_in_semantic_group",
    "n_leading_edge_genes",
    "top10_leading_edge_genes"
)

cat(
    "\n============================================================\n"
)

cat(
    "Positive: redundancy-reduced representative GO BP terms\n"
)

cat(
    "============================================================\n"
)

print(
    pos_summary[
        ,
        ..show_cols
    ]
)

cat(
    "\n============================================================\n"
)

cat(
    "Negative: redundancy-reduced representative GO BP terms\n"
)

cat(
    "============================================================\n"
)

print(
    neg_summary[
        ,
        ..show_cols
    ]
)


# ============================================================
# 13. Redundancy reduction statistics
# ============================================================

cat(
    "\nRedundancy reduction summary:\n"
)

cat(
    "Positive significant terms : ",
    nrow(pos),
    "\n",
    sep = ""
)

cat(
    "Positive representatives   : ",
    nrow(pos_red$representatives),
    "\n",
    sep = ""
)

cat(
    "Negative significant terms : ",
    nrow(neg),
    "\n",
    sep = ""
)

cat(
    "Negative representatives   : ",
    nrow(neg_red$representatives),
    "\n",
    sep = ""
)

cat(
    "Semantic similarity cutoff : ",
    SEMANTIC_CUTOFF,
    "\n",
    sep = ""
)

cat(
    "Semantic similarity method : ",
    SEMANTIC_METHOD,
    "\n",
    sep = ""
)


# ============================================================
# 14. Save complete bundle
# ============================================================

result <- list(

    parameters = list(
        FDR_cutoff =
            FDR_CUTOFF,

        semantic_cutoff =
            SEMANTIC_CUTOFF,

        semantic_method =
            SEMANTIC_METHOD,

        ontology =
            "BP",

        top_representative =
            TOP_REPRESENTATIVE
    ),

    positive = list(
        grouped_terms =
            pos_red$grouped,

        all_representatives =
            pos_red$representatives,

        top_representatives =
            pos_summary,

        leading_edge_long =
            pos_leading,

        semantic_similarity =
            pos_red$similarity
    ),

    negative = list(
        grouped_terms =
            neg_red$grouped,

        all_representatives =
            neg_red$representatives,

        top_representatives =
            neg_summary,

        leading_edge_long =
            neg_leading,

        semantic_similarity =
            neg_red$similarity
    )
)

saveRDS(
    result,
    OUT_RDS
)


# ============================================================
# 15. Final files
# ============================================================

cat(
    "\nSaved files:\n",
    OUT_POS_GROUPED, "\n",
    OUT_NEG_GROUPED, "\n",
    OUT_POS_REP_ALL, "\n",
    OUT_NEG_REP_ALL, "\n",
    OUT_POS_TOP, "\n",
    OUT_NEG_TOP, "\n",
    OUT_POS_LEADING_LONG, "\n",
    OUT_NEG_LEADING_LONG, "\n",
    OUT_POS_SIM, "\n",
    OUT_NEG_SIM, "\n",
    OUT_RDS, "\n",
    sep = ""
)
