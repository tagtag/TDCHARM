# ============================================================
# 14_core12_mouse_GO_BP_preranked_GSEA.R
# GSE303006 / GENCODE M23
#
# Preranked GSEA of the degree-corrected all-gene ranking
# for core #12 against Mouse Gene Ontology Biological Process.
#
# Input:
#   GSE303006_M23_core12_all16239_gene_rank_degree_corrected.rnk
#
# Ranking:
#   high positive score  -> positive side of core #12 E-P axis
#   strong negative score -> negative side of core #12 E-P axis
#
# Interpretation:
#   NES > 0 : GO BP enriched toward positive-ranked genes
#   NES < 0 : GO BP enriched toward negative-ranked genes
#
# Uses:
#   clusterProfiler::gseGO()
#   org.Mm.eg.db
#
# The script:
#   1) reads and validates the 16,239-gene ranked list
#   2) checks SYMBOL mapping to org.Mm.eg.db
#   3) runs preranked GO-BP GSEA
#   4) saves all GO terms
#   5) saves positive-NES and negative-NES tables separately
#   6) saves top significant terms for each side
# ============================================================


# ============================================================
# 0. Settings
# ============================================================

RNK_FILE <-
    "GSE303006_M23_core12_all16239_gene_rank_degree_corrected.rnk"

TOP_N <- 30L

FDR_CUTOFF <- 0.05

MIN_GS_SIZE <- 10L
MAX_GS_SIZE <- 500L

SEED <- 12345L


# Output files
OUT_MAPPING_QC <-
    "GSE303006_M23_core12_GO_BP_gene_mapping_QC.tsv"

OUT_ALL <-
    "GSE303006_M23_core12_GO_BP_GSEA_all.tsv"

OUT_POS <-
    "GSE303006_M23_core12_GO_BP_GSEA_positive_NES.tsv"

OUT_NEG <-
    "GSE303006_M23_core12_GO_BP_GSEA_negative_NES.tsv"

OUT_POS_SIG <-
    "GSE303006_M23_core12_GO_BP_GSEA_positive_FDR005.tsv"

OUT_NEG_SIG <-
    "GSE303006_M23_core12_GO_BP_GSEA_negative_FDR005.tsv"

OUT_POS_TOP <-
    "GSE303006_M23_core12_GO_BP_GSEA_positive_top30.tsv"

OUT_NEG_TOP <-
    "GSE303006_M23_core12_GO_BP_GSEA_negative_top30.tsv"

OUT_RDS <-
    "GSE303006_M23_core12_GO_BP_GSEA.rds"


# ============================================================
# 1. Package checks
# ============================================================

needed_packages <- c(
    "data.table",
    "clusterProfiler",
    "org.Mm.eg.db",
    "AnnotationDbi"
)

missing_packages <- needed_packages[
    !vapply(
        needed_packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]

if (length(missing_packages) > 0L) {

    stop(
        paste0(
            "Missing package(s): ",
            paste(
                missing_packages,
                collapse = ", "
            ),
            "\n\nInstall Bioconductor packages with:\n",
            "if (!requireNamespace(\"BiocManager\", quietly=TRUE)) ",
            "install.packages(\"BiocManager\")\n",
            "BiocManager::install(c(",
            paste(
                sprintf(
                    "\"%s\"",
                    setdiff(
                        missing_packages,
                        "data.table"
                    )
                ),
                collapse = ", "
            ),
            "))\n",
            if ("data.table" %in% missing_packages) {
                "install.packages(\"data.table\")\n"
            } else {
                ""
            }
        )
    )
}

library(data.table)


# ============================================================
# 2. Read preranked gene list
# ============================================================

if (!file.exists(RNK_FILE)) {
    stop(
        "Missing ranked-list file: ",
        RNK_FILE
    )
}

rank_dt <- fread(
    RNK_FILE,
    header = FALSE,
    col.names = c(
        "gene_name",
        "score"
    )
)

if (ncol(rank_dt) != 2L) {
    stop(
        "Expected exactly two columns in .rnk file."
    )
}

rank_dt[
    ,
    gene_name :=
        as.character(
            gene_name
        )
]

rank_dt[
    ,
    score :=
        as.numeric(
            score
        )
]

if (anyNA(rank_dt$gene_name) ||
    anyNA(rank_dt$score)) {
    stop(
        "NA found in gene name or ranking score."
    )
}

if (any(!is.finite(rank_dt$score))) {
    stop(
        "Non-finite ranking score found."
    )
}

cat(
    "Input genes :",
    nrow(rank_dt),
    "\n"
)


# ============================================================
# 3. Handle duplicated gene symbols defensively
#
# Expected here: no duplicated symbols.
# If duplicates exist, retain the row with the largest |score|.
# ============================================================

n_dup <- sum(
    duplicated(
        rank_dt$gene_name
    )
)

cat(
    "Duplicated gene symbols :",
    n_dup,
    "\n"
)

if (n_dup > 0L) {

    rank_dt[
        ,
        abs_score :=
            abs(
                score
            )
    ]

    setorder(
        rank_dt,
        gene_name,
        -abs_score
    )

    rank_dt <- rank_dt[
        ,
        .SD[1],
        by = gene_name
    ]

    rank_dt[
        ,
        abs_score := NULL
    ]
}


# ============================================================
# 4. Create named decreasing geneList
# ============================================================

setorder(
    rank_dt,
    -score,
    gene_name
)

geneList <- rank_dt$score

names(geneList) <- rank_dt$gene_name

if (!all(
    diff(geneList) <= 0
)) {
    stop(
        "geneList is not sorted decreasingly."
    )
}

if (anyDuplicated(
    names(geneList)
)) {
    stop(
        "Duplicated names remain in geneList."
    )
}

cat(
    "Genes after duplicate handling :",
    length(geneList),
    "\n"
)

cat(
    "Positive scores :",
    sum(geneList > 0),
    "\n"
)

cat(
    "Negative scores :",
    sum(geneList < 0),
    "\n"
)

cat(
    "Score range :",
    min(geneList),
    "to",
    max(geneList),
    "\n"
)


# ============================================================
# 5. Mouse SYMBOL annotation QC
#
# Check whether each ranked symbol exists in org.Mm.eg.db.
# This is annotation QC only; gseGO() will internally use
# the OrgDb with keyType="SYMBOL".
# ============================================================

mouse_symbols <- AnnotationDbi::keys(
    org.Mm.eg.db::org.Mm.eg.db,
    keytype = "SYMBOL"
)

mapping_qc <- data.table(
    gene_name =
        names(geneList),
    score =
        as.numeric(
            geneList
        ),
    found_in_orgMm =
        names(geneList) %in%
            mouse_symbols
)

cat(
    "\nGenes found in org.Mm.eg.db :",
    sum(mapping_qc$found_in_orgMm),
    "/",
    nrow(mapping_qc),
    "\n"
)

cat(
    "Mapping fraction :",
    mean(mapping_qc$found_in_orgMm),
    "\n"
)

fwrite(
    mapping_qc,
    OUT_MAPPING_QC,
    sep = "\t"
)


# ============================================================
# 6. Run Mouse GO Biological Process preranked GSEA
#
# pvalueCutoff = 1 is intentional:
# retain every tested term in the result and apply FDR filtering
# ourselves afterwards.
#
# method = "multilevel" is the current clusterProfiler default
# and supports preranked GSEA.
# ============================================================

set.seed(
    SEED
)

if (
    exists(
        "gsea_bp",
        inherits = TRUE
    ) &&
    inherits(
        get(
            "gsea_bp",
            inherits = TRUE
        ),
        "gseaResult"
    )
) {

    cat(
        "\nReusing existing gsea_bp object from the current R session.\n"
    )

    gsea_bp <- get(
        "gsea_bp",
        inherits = TRUE
    )

} else {

    cat(
        "\nRunning gseGO for Mouse GO Biological Process...\n"
    )

    gsea_bp <- clusterProfiler::gseGO(
        geneList =
            geneList,

        ont =
            "BP",

        OrgDb =
            org.Mm.eg.db::org.Mm.eg.db,

        keyType =
            "SYMBOL",

        exponent =
            1,

        minGSSize =
            MIN_GS_SIZE,

        maxGSSize =
            MAX_GS_SIZE,

        pvalueCutoff =
            1,

        pAdjustMethod =
            "BH",

        verbose =
            TRUE
    )
}

# ============================================================
# 7. Convert result to data.table
# ============================================================

res <- as.data.table(
    as.data.frame(
        gsea_bp
    )
)

cat(
    "\nGO BP terms tested/returned :",
    nrow(res),
    "\n"
)

if (nrow(res) == 0L) {
    stop(
        "gseGO returned zero GO BP terms."
    )
}

required_res_cols <- c(
    "ID",
    "Description",
    "setSize",
    "enrichmentScore",
    "NES",
    "pvalue",
    "p.adjust",
    "core_enrichment"
)

missing_res_cols <- setdiff(
    required_res_cols,
    names(res)
)

if (length(missing_res_cols) > 0L) {
    stop(
        "Unexpected gseGO result format. Missing: ",
        paste(
            missing_res_cols,
            collapse = ", "
        )
    )
}


# ============================================================
# 8. Add useful summary columns
# ============================================================

res[
    ,
    direction :=
        fifelse(
            NES > 0,
            "positive",
            fifelse(
                NES < 0,
                "negative",
                "zero"
            )
        )
]

res[
    ,
    significant_FDR005 :=
        p.adjust <=
        FDR_CUTOFF
]

# Leading-edge/core enrichment count
res[
    ,
    n_core_enrichment :=
        fifelse(
            is.na(
                core_enrichment
            ) |
            core_enrichment == "",
            0L,
            lengths(
                strsplit(
                    core_enrichment,
                    "/",
                    fixed = TRUE
                )
            )
        )
]

# Make primary ordering explicit:
# lowest FDR first, then strongest |NES|.
res[
    ,
    abs_NES :=
        abs(
            NES
        )
]

setorder(
    res,
    p.adjust,
    -abs_NES
)

fwrite(
    res,
    OUT_ALL,
    sep = "\t"
)


# ============================================================
# 9. Positive and negative NES tables
# ============================================================

pos <- res[
    NES > 0
][
    order(
        p.adjust,
        -NES
    )
]

neg <- res[
    NES < 0
][
    order(
        p.adjust,
        NES
    )
]

pos_sig <- pos[
    p.adjust <=
        FDR_CUTOFF
]

neg_sig <- neg[
    p.adjust <=
        FDR_CUTOFF
]

fwrite(
    pos,
    OUT_POS,
    sep = "\t"
)

fwrite(
    neg,
    OUT_NEG,
    sep = "\t"
)

fwrite(
    pos_sig,
    OUT_POS_SIG,
    sep = "\t"
)

fwrite(
    neg_sig,
    OUT_NEG_SIG,
    sep = "\t"
)


# ============================================================
# 10. Top positive / negative tables
#
# For the display tables:
#   primarily sort by FDR, then by NES magnitude.
#
# We keep up to TOP_N terms from each sign.
# ============================================================

pos_top <- pos[
    seq_len(
        min(
            TOP_N,
            .N
        )
    )
]

neg_top <- neg[
    seq_len(
        min(
            TOP_N,
            .N
        )
    )
]

fwrite(
    pos_top,
    OUT_POS_TOP,
    sep = "\t"
)

fwrite(
    neg_top,
    OUT_NEG_TOP,
    sep = "\t"
)


# ============================================================
# 11. Console tables
# ============================================================

show_cols <- c(
    "ID",
    "Description",
    "setSize",
    "NES",
    "pvalue",
    "p.adjust",
    "n_core_enrichment"
)

cat(
    "\n============================================================\n"
)

cat(
    "Positive NES: top GO Biological Process terms\n"
)

cat(
    "============================================================\n"
)

print(
    pos_top[
        ,
        ..show_cols
    ]
)

cat(
    "\n============================================================\n"
)

cat(
    "Negative NES: top GO Biological Process terms\n"
)

cat(
    "============================================================\n"
)

print(
    neg_top[
        ,
        ..show_cols
    ]
)


# ============================================================
# 12. Significant-term summary
# ============================================================

cat(
    "\nSignificant GO BP terms at FDR <=",
    FDR_CUTOFF,
    ":\n"
)

cat(
    "Positive NES :",
    nrow(pos_sig),
    "\n"
)

cat(
    "Negative NES :",
    nrow(neg_sig),
    "\n"
)

if (nrow(pos_sig) > 0L) {

    cat(
        "\nTop significant positive NES terms:\n"
    )

    print(
        pos_sig[
            1:min(
                20L,
                .N
            ),
            ..show_cols
        ]
    )
}

if (nrow(neg_sig) > 0L) {

    cat(
        "\nTop significant negative NES terms:\n"
    )

    print(
        neg_sig[
            1:min(
                20L,
                .N
            ),
            ..show_cols
        ]
    )
}


# ============================================================
# 13. Save complete R object
# ============================================================

result <- list(

    input_file =
        RNK_FILE,

    parameters = list(
        ontology =
            "BP",

        species =
            "Mus musculus",

        keyType =
            "SYMBOL",

        exponent =
            1,

        minGSSize =
            MIN_GS_SIZE,

        maxGSSize =
            MAX_GS_SIZE,

        FDR_cutoff =
            FDR_CUTOFF,

        method =
            "clusterProfiler/fgsea default",

        seed =
            SEED
    ),

    mapping_QC =
        mapping_qc,

    gsea_object =
        gsea_bp,

    all_results =
        res,

    positive_results =
        pos,

    negative_results =
        neg,

    positive_significant =
        pos_sig,

    negative_significant =
        neg_sig,

    positive_top =
        pos_top,

    negative_top =
        neg_top
)

saveRDS(
    result,
    OUT_RDS
)


# ============================================================
# 14. Final summary
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Mouse GO BP preranked GSEA completed\n"
)

cat(
    "============================================================\n"
)

cat(
    "Ranked genes :",
    length(geneList),
    "\n"
)

cat(
    "Mapped SYMBOLs :",
    sum(mapping_qc$found_in_orgMm),
    "\n"
)

cat(
    "Returned GO BP terms :",
    nrow(res),
    "\n"
)

cat(
    "FDR <= 0.05 positive NES :",
    nrow(pos_sig),
    "\n"
)

cat(
    "FDR <= 0.05 negative NES :",
    nrow(neg_sig),
    "\n"
)

cat(
    "\nSaved files:\n",
    OUT_MAPPING_QC, "\n",
    OUT_ALL, "\n",
    OUT_POS, "\n",
    OUT_NEG, "\n",
    OUT_POS_SIG, "\n",
    OUT_NEG_SIG, "\n",
    OUT_POS_TOP, "\n",
    OUT_NEG_TOP, "\n",
    OUT_RDS, "\n",
    sep = ""
)
