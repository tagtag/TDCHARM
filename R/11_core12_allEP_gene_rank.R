# ============================================================
# 13_core12_allEP_gene_rank_degree_corrected.R
# GSE303006 / GENCODE M23
#
# Build a DEGREE-CORRECTED ranked list of all 16,239 target genes
# from ALL 730,969 E-P pairs for core rank 12.
#
# Core-specific E-P score:
#
#   s_p = CORE[a,b,k,r] * U_EP[p,a]
#
# where (a,b,k,r) are taken from core rank 12.
#
# Gene-level degree-corrected score:
#
#   S_g = (1 / sqrt(d_g)) * sum_{p:g(p)=g} s_p
#
# This is the projection of the E-P score vector onto the
# normalized incidence vector of each gene. If a gene-level
# signal q_g were replicated over d_g E-P pairs as
# q_g / sqrt(d_g), this aggregation recovers q_g exactly.
#
# We also save diagnostic quantities:
#   raw sum
#   mean
#   positive contribution
#   negative contribution
#   absolute contribution
#   max |E-P score|
#   number of E-P pairs / enhancers
#   Stage-1 RNA score at the original K used by core #12
#
# Primary GSEA rank:
#   gene_score_degree_corrected
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 0. Settings / files
# ------------------------------------------------------------

CORE_RANK <- 12L

CORE_FILE <-
    "GSE303006_M23_stage2_Tucker_core.rds"

TOP_CORE_FILE <-
    "GSE303006_M23_stage2_top100_core.tsv"

EP_FACTOR_FILE <-
    "GSE303006_M23_stage2_Tucker_EP_factor.rds"

MAPPING_FILE <-
    "GSE303006_M23_stage2_EP_mapping_degrees.rds"

RNA_STAGE1_FILE <-
    "GSE303006_M23_stage1_K20_RNA_signaligned.rds"

OUT_GENE_TABLE <-
    "GSE303006_M23_core12_all16239_gene_rank_degree_corrected.tsv"

OUT_GSEA_RNK <-
    "GSE303006_M23_core12_all16239_gene_rank_degree_corrected.rnk"

OUT_POSITIVE <-
    "GSE303006_M23_core12_gene_rank_positive.tsv"

OUT_NEGATIVE <-
    "GSE303006_M23_core12_gene_rank_negative.tsv"

OUT_EP_SCORE <-
    "GSE303006_M23_core12_all730969_EP_scores.rds"

OUT_RDS <-
    "GSE303006_M23_core12_gene_rank_degree_corrected.rds"


# ------------------------------------------------------------
# 1. Input checks
# ------------------------------------------------------------

need <- c(
    CORE_FILE,
    TOP_CORE_FILE,
    EP_FACTOR_FILE,
    MAPPING_FILE,
    RNA_STAGE1_FILE
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
# 2. Load core #12 information
# ------------------------------------------------------------

top_core <- fread(
    TOP_CORE_FILE
)

core_info <- top_core[
    core_rank == CORE_RANK
]

if (nrow(core_info) != 1L) {
    stop(
        "Expected exactly one row for core rank ",
        CORE_RANK,
        "; found ",
        nrow(core_info)
    )
}

a <- as.integer(
    core_info$a_EP
)

b <- as.integer(
    core_info$b_CELL
)

k <- as.integer(
    core_info$k
)

r <- as.integer(
    core_info$r_MODAL
)

core_value_from_table <- as.numeric(
    core_info$core_value
)

CORE <- readRDS(
    CORE_FILE
)

if (length(dim(CORE)) != 4L) {
    stop(
        "CORE is not a four-way array."
    )
}

if (
    a < 1L || a > dim(CORE)[1] ||
    b < 1L || b > dim(CORE)[2] ||
    k < 1L || k > dim(CORE)[3] ||
    r < 1L || r > dim(CORE)[4]
) {
    stop(
        "Core indices are out of range."
    )
}

core_value <- as.numeric(
    CORE[
        a,
        b,
        k,
        r
    ]
)

if (
    abs(
        core_value -
        core_value_from_table
    ) > 1e-12
) {
    warning(
        "core_value in table and CORE differ slightly."
    )
}

cat(
    "Core rank       : ",
    CORE_RANK,
    "\n",
    sep = ""
)

cat(
    "(a,b,k,r)       : (",
    a, ", ",
    b, ", ",
    k, ", ",
    r, ")\n",
    sep = ""
)

cat(
    "CORE value      : ",
    format(
        core_value,
        digits = 12
    ),
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 3. Load E-P factor and mapping
# ------------------------------------------------------------

U_EP <- readRDS(
    EP_FACTOR_FILE
)

map <- as.data.table(
    readRDS(
        MAPPING_FILE
    )
)

required_map <- c(
    "stage2_pair_id",
    "enhancer_id",
    "gene_name",
    "binpair_id",
    "d_g"
)

if (!all(required_map %in% names(map))) {
    stop(
        "Mapping lacks required columns: ",
        paste(
            setdiff(
                required_map,
                names(map)
            ),
            collapse = ", "
        )
    )
}

if (nrow(U_EP) != nrow(map)) {
    stop(
        "U_EP rows and mapping rows differ: ",
        nrow(U_EP),
        " vs ",
        nrow(map)
    )
}

if (a > ncol(U_EP)) {
    stop(
        "Requested E-P component exceeds U_EP columns."
    )
}

if (!is.null(rownames(U_EP))) {

    if (!identical(
        rownames(U_EP),
        as.character(
            map$stage2_pair_id
        )
    )) {
        stop(
            "U_EP row order does not match stage2_pair_id."
        )
    }
}

cat(
    "E-P pairs       : ",
    nrow(map),
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 4. Score ALL 730,969 E-P pairs for core #12
#
# s_p = core_value * U_EP[p,a]
#
# This is the signed E-P-side amplitude associated with
# this individual core element.
# ------------------------------------------------------------

ep_loading <- U_EP[
    ,
    a
]

ep_score <- core_value *
    ep_loading

if (any(!is.finite(ep_score))) {
    stop(
        "Non-finite E-P scores found."
    )
}

ep_score_dt <- data.table(
    stage2_pair_id =
        map$stage2_pair_id,
    enhancer_id =
        map$enhancer_id,
    gene_name =
        map$gene_name,
    binpair_id =
        map$binpair_id,
    d_g =
        map$d_g,
    EP_loading =
        ep_loading,
    core_EP_score =
        ep_score
)

# Save compact full E-P score object as RDS.
# TSV is intentionally not written by default because this is 730,969 rows.
saveRDS(
    ep_score_dt,
    OUT_EP_SCORE
)

cat(
    "\nAll-E-P score summary:\n"
)

print(
    summary(
        ep_score
    )
)


# ------------------------------------------------------------
# 5. Verify d_g from the actual final E-P mapping
#
# The stored d_g should exactly equal the number of final
# E-P pairs mapped to each gene.
# ------------------------------------------------------------

degree_check <- map[
    ,
    .(
        d_g_recount = .N,
        d_g_stored =
            unique(
                d_g
            )
    ),
    by = gene_name
]

if (
    any(
        lengths(
            map[
                ,
                .(
                    vals = list(
                        unique(
                            d_g
                        )
                    )
                ),
                by = gene_name
            ]$vals
        ) != 1L
    )
) {
    stop(
        "d_g is not constant within gene."
    )
}

if (
    any(
        degree_check$d_g_recount !=
            degree_check$d_g_stored
    )
) {
    stop(
        "Stored d_g does not match E-P counts."
    )
}

cat(
    "\nd_g verification passed.\n"
)


# ------------------------------------------------------------
# 6. Gene-level aggregation from ALL E-P pairs
#
# Primary score:
#
#   S_g = sum(s_p) / sqrt(d_g)
#
# Additional diagnostics:
#
#   positive_score_degree_corrected
#     = sum(max(s_p,0)) / sqrt(d_g)
#
#   negative_score_degree_corrected
#     = sum(min(s_p,0)) / sqrt(d_g)
#
#   absolute_score_degree_corrected
#     = sum(|s_p|) / sqrt(d_g)
#
# Note:
# positive + negative = net score.
# ------------------------------------------------------------

gene_rank <- ep_score_dt[
    ,
    .(
        d_g = first(d_g),

        n_EP = .N,

        n_unique_enhancers =
            uniqueN(
                enhancer_id
            ),

        raw_sum_EP_score =
            sum(
                core_EP_score
            ),

        mean_EP_score =
            mean(
                core_EP_score
            ),

        max_EP_score =
            max(
                core_EP_score
            ),

        min_EP_score =
            min(
                core_EP_score
            ),

        max_abs_EP_score =
            max(
                abs(
                    core_EP_score
                )
            ),

        sum_abs_EP_score =
            sum(
                abs(
                    core_EP_score
                )
            ),

        positive_raw_sum =
            sum(
                pmax(
                    core_EP_score,
                    0
                )
            ),

        negative_raw_sum =
            sum(
                pmin(
                    core_EP_score,
                    0
                )
            ),

        n_positive_EP =
            sum(
                core_EP_score > 0
            ),

        n_negative_EP =
            sum(
                core_EP_score < 0
            )
    ),
    by = gene_name
]

if (
    any(
        gene_rank$n_EP !=
            gene_rank$d_g
    )
) {
    stop(
        "Gene E-P counts and d_g disagree."
    )
}

gene_rank[
    ,
    gene_score_degree_corrected :=
        raw_sum_EP_score /
        sqrt(
            d_g
        )
]

gene_rank[
    ,
    positive_score_degree_corrected :=
        positive_raw_sum /
        sqrt(
            d_g
        )
]

gene_rank[
    ,
    negative_score_degree_corrected :=
        negative_raw_sum /
        sqrt(
            d_g
        )
]

gene_rank[
    ,
    absolute_score_degree_corrected :=
        sum_abs_EP_score /
        sqrt(
            d_g
        )
]


# ------------------------------------------------------------
# 7. Add Stage-1 RNA score at core-specific K
#
# This is NOT used for the primary gene ranking.
# It is provided as biological/latent-component annotation.
# ------------------------------------------------------------

rna_obj <- readRDS(
    RNA_STAGE1_FILE
)

if (!"feature_scores" %in% names(rna_obj)) {
    stop(
        "RNA Stage-1 object lacks feature_scores."
    )
}

RNA_B <- rna_obj$feature_scores

if (k > ncol(RNA_B)) {
    stop(
        "Core k exceeds RNA feature-score columns."
    )
}

rna_k_score <- RNA_B[
    ,
    k
]

names(rna_k_score) <- rownames(
    RNA_B
)

gene_rank[
    ,
    RNA_K_score :=
        unname(
            rna_k_score[
                gene_name
            ]
        )
]

gene_rank[
    ,
    RNA_K_abs_score :=
        abs(
            RNA_K_score
        )
]

if (anyNA(
    gene_rank$RNA_K_score
)) {

    warning(
        sum(
            is.na(
                gene_rank$RNA_K_score
            )
        ),
        " genes have no RNA K score."
    )
}


# ------------------------------------------------------------
# 8. Ranking
#
# Signed ranking:
#   large positive -> top
#   large negative -> bottom
#
# Absolute ranking is also saved for diagnostics.
# ------------------------------------------------------------

setorder(
    gene_rank,
    -gene_score_degree_corrected,
    gene_name
)

gene_rank[
    ,
    signed_rank :=
        seq_len(.N)
]

gene_rank[
    ,
    abs_rank :=
        frank(
            -abs(
                gene_score_degree_corrected
            ),
            ties.method = "first"
        )
]

gene_rank[
    ,
    positive_rank :=
        fifelse(
            gene_score_degree_corrected > 0,
            frank(
                -gene_score_degree_corrected,
                ties.method = "first",
                na.last = "keep"
            ),
            NA_integer_
        )
]

gene_rank[
    ,
    negative_rank :=
        fifelse(
            gene_score_degree_corrected < 0,
            frank(
                gene_score_degree_corrected,
                ties.method = "first",
                na.last = "keep"
            ),
            NA_integer_
        )
]


# ------------------------------------------------------------
# 9. Basic checks and summaries
# ------------------------------------------------------------

n_genes <- nrow(
    gene_rank
)

cat(
    "\nUnique genes in final E-P set : ",
    n_genes,
    "\n",
    sep = ""
)

if (n_genes != 16239L) {
    warning(
        "Expected 16,239 genes but obtained ",
        n_genes
    )
}

cat(
    "Positive gene scores          : ",
    sum(
        gene_rank$gene_score_degree_corrected > 0
    ),
    "\n",
    sep = ""
)

cat(
    "Negative gene scores          : ",
    sum(
        gene_rank$gene_score_degree_corrected < 0
    ),
    "\n",
    sep = ""
)

cat(
    "Zero gene scores              : ",
    sum(
        gene_rank$gene_score_degree_corrected == 0
    ),
    "\n",
    sep = ""
)

cat(
    "\nDegree-corrected gene-score summary:\n"
)

print(
    summary(
        gene_rank$gene_score_degree_corrected
    )
)


# ------------------------------------------------------------
# 10. Save full table and GSEA .rnk
# ------------------------------------------------------------

fwrite(
    gene_rank,
    OUT_GENE_TABLE,
    sep = "\t"
)

# GSEA preranked format: gene symbol + signed score.
# No header for broad compatibility.
fwrite(
    gene_rank[
        ,
        .(
            gene_name,
            gene_score_degree_corrected
        )
    ],
    OUT_GSEA_RNK,
    sep = "\t",
    col.names = FALSE
)


# ------------------------------------------------------------
# 11. Positive / negative ranked tables
# ------------------------------------------------------------

positive_genes <- gene_rank[
    gene_score_degree_corrected > 0
][
    order(
        -gene_score_degree_corrected
    )
]

negative_genes <- gene_rank[
    gene_score_degree_corrected < 0
][
    order(
        gene_score_degree_corrected
    )
]

fwrite(
    positive_genes,
    OUT_POSITIVE,
    sep = "\t"
)

fwrite(
    negative_genes,
    OUT_NEGATIVE,
    sep = "\t"
)


# ------------------------------------------------------------
# 12. Console: top/bottom 30
# ------------------------------------------------------------

show_cols <- c(
    "signed_rank",
    "gene_name",
    "d_g",
    "n_unique_enhancers",
    "gene_score_degree_corrected",
    "positive_score_degree_corrected",
    "negative_score_degree_corrected",
    "max_abs_EP_score",
    "RNA_K_score"
)

cat(
    "\nTop 30 positive-ranked genes:\n"
)

print(
    gene_rank[
        1:min(
            30L,
            .N
        ),
        ..show_cols
    ]
)

cat(
    "\nTop 30 negative-ranked genes:\n"
)

print(
    gene_rank[
        (.N - min(30L, .N) + 1L):.N,
        ..show_cols
    ][
        order(
            gene_score_degree_corrected
        )
    ]
)


# ------------------------------------------------------------
# 13. Degree-bias diagnostics
#
# The point of degree correction is to prevent genes with many
# E-P pairs from dominating merely because d_g is large.
# We therefore report Spearman correlations with d_g.
# ------------------------------------------------------------

cor_raw_degree <- suppressWarnings(
    cor(
        gene_rank$d_g,
        abs(
            gene_rank$raw_sum_EP_score
        ),
        method = "spearman"
    )
)

cor_corrected_degree <- suppressWarnings(
    cor(
        gene_rank$d_g,
        abs(
            gene_rank$gene_score_degree_corrected
        ),
        method = "spearman"
    )
)

cat(
    "\nDegree-bias diagnostic (Spearman rho):\n"
)

cat(
    "d_g vs |raw sum|             : ",
    cor_raw_degree,
    "\n",
    sep = ""
)

cat(
    "d_g vs |degree-corrected|    : ",
    cor_corrected_degree,
    "\n",
    sep = ""
)


# ------------------------------------------------------------
# 14. Save RDS bundle
# ------------------------------------------------------------

result <- list(

    core_rank =
        CORE_RANK,

    core_indices = c(
        a_EP = a,
        b_CELL = b,
        k = k,
        r_MODAL = r
    ),

    core_value =
        core_value,

    gene_score_definition =
        "sum(core_value * U_EP[p,a] over E-P pairs targeting gene g) / sqrt(d_g)",

    n_EP =
        nrow(map),

    n_genes =
        n_genes,

    gene_rank =
        gene_rank,

    degree_bias = list(
        spearman_dg_vs_abs_raw_sum =
            cor_raw_degree,

        spearman_dg_vs_abs_degree_corrected =
            cor_corrected_degree
    ),

    files = list(
        full_gene_table =
            OUT_GENE_TABLE,

        GSEA_rank =
            OUT_GSEA_RNK,

        positive_genes =
            OUT_POSITIVE,

        negative_genes =
            OUT_NEGATIVE,

        full_EP_scores =
            OUT_EP_SCORE
    )
)

saveRDS(
    result,
    OUT_RDS
)


# ------------------------------------------------------------
# 15. Final message
# ------------------------------------------------------------

cat(
    "\n============================================================\n"
)

cat(
    "Core #12 all-E-P -> gene ranking completed\n"
)

cat(
    "============================================================\n"
)

cat(
    "E-P pairs used : ",
    nrow(map),
    "\n",
    sep = ""
)

cat(
    "Genes ranked   : ",
    n_genes,
    "\n",
    sep = ""
)

cat(
    "Primary rank   : sum(E-P score) / sqrt(d_g)\n"
)

cat(
    "\nSaved files:\n",
    OUT_GENE_TABLE, "\n",
    OUT_GSEA_RNK, "\n",
    OUT_POSITIVE, "\n",
    OUT_NEGATIVE, "\n",
    OUT_EP_SCORE, "\n",
    OUT_RDS, "\n",
    sep = ""
)
