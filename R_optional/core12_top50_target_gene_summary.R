# ============================================================
# 12_core12_target_gene_summary.R
# GSE303006 / GENCODE M23
#
# Gene-level characterization of target genes from:
#   core rank 12
#   positive-loading top 50 E-P pairs
#   negative-loading top 50 E-P pairs
#
# Outputs:
#   - the 100 selected E-P pairs
#   - gene-level summary for positive and negative sides
#   - genes appearing on both sides
#   - Stage-1 RNA K2 feature score for each target gene
# ============================================================

library(data.table)

CORE_RANK <- 12L
TOP_N <- 50L

TOP_EP_FILE <-
    "GSE303006_M23_stage2_top100_core_top_EP.tsv.gz"

# Accept uncompressed file as fallback.
if (!file.exists(TOP_EP_FILE)) {
    alt <- sub("\\.gz$", "", TOP_EP_FILE)
    if (file.exists(alt)) {
        TOP_EP_FILE <- alt
    }
}

RNA_STAGE1_FILE <-
    "GSE303006_M23_stage1_K20_RNA_signaligned.rds"

TOP_CORE_FILE <-
    "GSE303006_M23_stage2_top100_core.tsv"

OUT_EP <-
    "GSE303006_M23_core12_top50pos_top50neg_EP.tsv"

OUT_GENE <-
    "GSE303006_M23_core12_target_gene_summary.tsv"

OUT_POS_GENE <-
    "GSE303006_M23_core12_positive_target_genes.tsv"

OUT_NEG_GENE <-
    "GSE303006_M23_core12_negative_target_genes.tsv"

OUT_OVERLAP <-
    "GSE303006_M23_core12_target_genes_both_sides.tsv"

OUT_RDS <-
    "GSE303006_M23_core12_target_gene_summary.rds"


# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

need <- c(
    TOP_EP_FILE,
    RNA_STAGE1_FILE,
    TOP_CORE_FILE
)

missing <- need[!file.exists(need)]

if (length(missing) > 0L) {
    stop(
        "Missing input file(s):\n",
        paste(missing, collapse = "\n")
    )
}


# ------------------------------------------------------------
# 2. Read selected E-P pairs and core information
# ------------------------------------------------------------

ep <- fread(TOP_EP_FILE)

required_cols <- c(
    "core_rank",
    "direction",
    "rank_within_direction",
    "enhancer_id",
    "gene_name",
    "signed_core_EP_score"
)

if (!all(required_cols %in% names(ep))) {
    stop(
        "Missing required columns in top E-P file: ",
        paste(
            setdiff(required_cols, names(ep)),
            collapse = ", "
        )
    )
}

core_info <- fread(TOP_CORE_FILE)[
    core_rank == CORE_RANK
]

if (nrow(core_info) != 1L) {
    stop(
        "Expected exactly one row for core rank ",
        CORE_RANK,
        ", found ",
        nrow(core_info)
    )
}

cat("\nCore #", CORE_RANK, ":\n", sep = "")
print(core_info)

# Restrict strictly to top 50 positive and top 50 negative.
ep12 <- ep[
    core_rank == CORE_RANK &
    direction %in% c("positive", "negative") &
    rank_within_direction <= TOP_N
]

setorder(
    ep12,
    direction,
    rank_within_direction
)

cat(
    "\nSelected E-P pairs : ",
    nrow(ep12),
    "\n",
    sep = ""
)

cat(
    "Positive E-P pairs : ",
    ep12[direction == "positive", .N],
    "\n",
    sep = ""
)

cat(
    "Negative E-P pairs : ",
    ep12[direction == "negative", .N],
    "\n",
    sep = ""
)

if (ep12[direction == "positive", .N] != TOP_N ||
    ep12[direction == "negative", .N] != TOP_N) {
    warning(
        "Expected 50 positive and 50 negative E-P pairs."
    )
}

fwrite(
    ep12,
    OUT_EP,
    sep = "\t"
)


# ------------------------------------------------------------
# 3. Load Stage-1 RNA feature scores
#
# Core #12 should have a definite original K, expected K2
# from the preceding analysis. Read k from top_core rather
# than hard-coding it.
# ------------------------------------------------------------

k_core <- as.integer(
    core_info$k
)

cat(
    "\nOriginal Stage-1 K used by core #",
    CORE_RANK,
    " : K",
    k_core,
    "\n",
    sep = ""
)

rna_obj <- readRDS(
    RNA_STAGE1_FILE
)

if (!"feature_scores" %in% names(rna_obj)) {
    stop(
        "RNA Stage-1 object has no feature_scores."
    )
}

RNA_B <- rna_obj$feature_scores

if (k_core < 1L || k_core > ncol(RNA_B)) {
    stop(
        "Invalid k in core information: ",
        k_core
    )
}

rna_k_score <- RNA_B[
    ,
    k_core
]

names(rna_k_score) <- rownames(
    RNA_B
)


# ------------------------------------------------------------
# 4. Gene-level aggregation separately by direction
#
# We keep several summaries:
#
# n_EP:
#   number of selected E-P pairs pointing to the gene
#
# n_unique_enhancers:
#   number of distinct selected enhancers targeting the gene
#
# sum_signed_EP_score:
#   total signed core-associated E-P score
#
# mean_signed_EP_score:
#   average signed E-P score
#
# max_abs_EP_score:
#   strongest single selected E-P for the gene
#
# sum_abs_EP_score:
#   total magnitude across selected E-P pairs
# ------------------------------------------------------------

gene_dir <- ep12[
    ,
    .(
        n_EP = .N,

        n_unique_enhancers =
            uniqueN(
                enhancer_id
            ),

        sum_signed_EP_score =
            sum(
                signed_core_EP_score
            ),

        mean_signed_EP_score =
            mean(
                signed_core_EP_score
            ),

        max_abs_EP_score =
            max(
                abs(
                    signed_core_EP_score
                )
            ),

        sum_abs_EP_score =
            sum(
                abs(
                    signed_core_EP_score
                )
            ),

        strongest_EP_stage2_pair_id =
            stage2_pair_id[
                which.max(
                    abs(
                        signed_core_EP_score
                    )
                )
            ],

        strongest_enhancer =
            enhancer_id[
                which.max(
                    abs(
                        signed_core_EP_score
                    )
                )
            ]
    ),
    by = .(
        direction,
        gene_name
    )
]

gene_dir[
    ,
    RNA_K_score :=
        unname(
            rna_k_score[
                gene_name
            ]
        )
]

gene_dir[
    ,
    RNA_K_abs_score :=
        abs(
            RNA_K_score
        )
]

# Rank genes separately within positive / negative side.
gene_dir[
    ,
    rank_by_sum_abs_EP :=
        frank(
            -sum_abs_EP_score,
            ties.method = "min"
        ),
    by = direction
]

gene_dir[
    ,
    rank_by_max_abs_EP :=
        frank(
            -max_abs_EP_score,
            ties.method = "min"
        ),
    by = direction
]

gene_dir[
    ,
    rank_by_RNA_K_abs :=
        frank(
            -RNA_K_abs_score,
            ties.method = "min",
            na.last = "keep"
        ),
    by = direction
]


# ------------------------------------------------------------
# 5. Mark genes appearing on both sides
# ------------------------------------------------------------

gene_presence <- gene_dir[
    ,
    .(
        sides =
            paste(
                sort(
                    unique(
                        direction
                    )
                ),
                collapse = "+"
            ),

        n_sides =
            uniqueN(
                direction
            )
    ),
    by = gene_name
]

gene_dir[
    gene_presence,
    on = "gene_name",
    `:=`(
        sides = i.sides,
        n_sides = i.n_sides
    )
]

gene_dir[
    ,
    appears_both_sides :=
        n_sides == 2L
]


# ------------------------------------------------------------
# 6. Write separate positive / negative lists
# ------------------------------------------------------------

pos_gene <- gene_dir[
    direction == "positive"
][
    order(
        -sum_abs_EP_score,
        -max_abs_EP_score
    )
]

neg_gene <- gene_dir[
    direction == "negative"
][
    order(
        -sum_abs_EP_score,
        -max_abs_EP_score
    )
]

fwrite(
    pos_gene,
    OUT_POS_GENE,
    sep = "\t"
)

fwrite(
    neg_gene,
    OUT_NEG_GENE,
    sep = "\t"
)

fwrite(
    gene_dir[
        order(
            direction,
            -sum_abs_EP_score
        )
    ],
    OUT_GENE,
    sep = "\t"
)


# ------------------------------------------------------------
# 7. Explicit positive/negative overlap table
# ------------------------------------------------------------

both_genes <- gene_presence[
    n_sides == 2L,
    gene_name
]

if (length(both_genes) > 0L) {

    overlap <- dcast(
        gene_dir[
            gene_name %in% both_genes
        ],
        gene_name + RNA_K_score ~ direction,
        value.var = c(
            "n_EP",
            "n_unique_enhancers",
            "sum_signed_EP_score",
            "sum_abs_EP_score",
            "max_abs_EP_score"
        )
    )

} else {

    overlap <- data.table(
        gene_name = character(),
        RNA_K_score = numeric()
    )
}

fwrite(
    overlap,
    OUT_OVERLAP,
    sep = "\t"
)


# ------------------------------------------------------------
# 8. Basic biological/statistical summaries
# ------------------------------------------------------------

n_pos_gene <- uniqueN(
    ep12[
        direction == "positive",
        gene_name
    ]
)

n_neg_gene <- uniqueN(
    ep12[
        direction == "negative",
        gene_name
    ]
)

n_both <- length(
    both_genes
)

cat(
    "\n============================================\n"
)

cat(
    "Core #",
    CORE_RANK,
    " target-gene summary\n",
    sep = ""
)

cat(
    "============================================\n"
)

cat(
    "Positive 50 E-P -> unique genes : ",
    n_pos_gene,
    "\n",
    sep = ""
)

cat(
    "Negative 50 E-P -> unique genes : ",
    n_neg_gene,
    "\n",
    sep = ""
)

cat(
    "Genes appearing on BOTH sides   : ",
    n_both,
    "\n",
    sep = ""
)

cat(
    "Union of target genes           : ",
    uniqueN(
        ep12$gene_name
    ),
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 9. Compare RNA K score between positive and negative genes
#
# This is descriptive because genes are selected through the
# same hierarchical analysis; it is not an independent test.
# ------------------------------------------------------------

pos_scores <- unique(
    pos_gene[
        !is.na(RNA_K_score),
        .(
            gene_name,
            RNA_K_score
        )
    ]
)$RNA_K_score

neg_scores <- unique(
    neg_gene[
        !is.na(RNA_K_score),
        .(
            gene_name,
            RNA_K_score
        )
    ]
)$RNA_K_score

cat(
    "\nRNA K",
    k_core,
    " scores of target genes:\n",
    sep = ""
)

cat(
    "Positive-side genes:\n"
)
print(
    summary(
        pos_scores
    )
)

cat(
    "Negative-side genes:\n"
)
print(
    summary(
        neg_scores
    )
)

if (length(pos_scores) > 0L &&
    length(neg_scores) > 0L) {

    wt <- wilcox.test(
        pos_scores,
        neg_scores,
        exact = FALSE
    )

    cat(
        "\nWilcoxon comparison of RNA K",
        k_core,
        " scores (descriptive):\n",
        sep = ""
    )

    print(wt)
}


# ------------------------------------------------------------
# 10. Top genes for console inspection
# ------------------------------------------------------------

cat(
    "\nTop positive-side target genes:\n"
)

print(
    pos_gene[
        1:min(
            20L,
            .N
        ),
        .(
            gene_name,
            n_EP,
            n_unique_enhancers,
            sum_signed_EP_score,
            sum_abs_EP_score,
            max_abs_EP_score,
            RNA_K_score
        )
    ]
)

cat(
    "\nTop negative-side target genes:\n"
)

print(
    neg_gene[
        1:min(
            20L,
            .N
        ),
        .(
            gene_name,
            n_EP,
            n_unique_enhancers,
            sum_signed_EP_score,
            sum_abs_EP_score,
            max_abs_EP_score,
            RNA_K_score
        )
    ]
)

if (n_both > 0L) {

    cat(
        "\nGenes appearing on both sides:\n"
    )

    print(
        overlap
    )
}


# ------------------------------------------------------------
# 11. Save RDS bundle
# ------------------------------------------------------------

result <- list(
    core_info =
        core_info,

    k_core =
        k_core,

    selected_EP =
        ep12,

    gene_summary =
        gene_dir,

    positive_genes =
        pos_gene,

    negative_genes =
        neg_gene,

    genes_both_sides =
        overlap
)

saveRDS(
    result,
    OUT_RDS
)


# ------------------------------------------------------------
# 12. Final file list
# ------------------------------------------------------------

cat(
    "\nSaved files:\n",
    OUT_EP, "\n",
    OUT_GENE, "\n",
    OUT_POS_GENE, "\n",
    OUT_NEG_GENE, "\n",
    OUT_OVERLAP, "\n",
    OUT_RDS, "\n",
    sep = ""
)
