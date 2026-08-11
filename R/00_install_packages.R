# ============================================================
# 00_install_packages.R
#
# Install/check R packages required by the GSE303006
# hierarchical tensor-decomposition workflow.
#
# External non-R requirements:
#   - MACS3 (only needed if peaks are regenerated)
#   - zcat / awk (used by preprocessing on Linux/Unix)
# ============================================================

cran <- c(
    "data.table",
    "Matrix",
    "RSpectra",
    "igraph"
)

bioc <- c(
    "GenomicRanges",
    "rtracklayer",
    "rhdf5",
    "clusterProfiler",
    "org.Mm.eg.db",
    "AnnotationDbi",
    "GOSemSim",
    "GO.db"
)

missing_cran <- cran[
    !vapply(cran, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_cran) > 0L) {
    install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

missing_bioc <- bioc[
    !vapply(bioc, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_bioc) > 0L) {
    BiocManager::install(
        missing_bioc,
        ask = FALSE,
        update = FALSE
    )
}

cat("\nPackage check completed.\n")
cat("R version:\n")
print(R.version.string)

cat("\nInstalled package versions:\n")
pkgs <- c(cran, bioc)
versions <- vapply(
    pkgs,
    function(p) as.character(utils::packageVersion(p)),
    character(1)
)
print(versions)
