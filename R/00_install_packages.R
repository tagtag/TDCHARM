# ============================================================
# 00_install_packages.R
#
# Install/check R packages required by the GSE303006
# hierarchical tensor-decomposition workflow.
#
# Usage:
#   Rscript R/00_install_packages.R                         # install missing
#   TDCHARM_CHECK_ONLY=1 Rscript R/00_install_packages.R    # check only
#
# Run it on a node with internet access (e.g. the login node), with
# the same R that the jobs use (e.g. after conda activate r_env).
#
# External non-R requirements:
#   - MACS3 (only needed if peaks are regenerated)
#   - zcat / awk (used by preprocessing on Linux/Unix)
# ============================================================

cran <- c(
    "data.table",
    "R.utils",      # needed by data.table::fread() to read .gz files
    "readxl",       # 23b: CHARM Supplementary Table 7 (.xlsx)
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

# A CRAN mirror must be set for non-interactive installation
# (Rscript); use the cloud mirror unless one is already configured.
repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
}

is_missing <- function(p) !requireNamespace(p, quietly = TRUE)

missing_cran <- cran[vapply(cran, is_missing, logical(1))]
missing_bioc <- bioc[vapply(bioc, is_missing, logical(1))]

cat("R: ", R.version.string, "\n", sep = "")
cat("Library paths:\n"); print(.libPaths())
cat("\nMissing CRAN packages        : ",
    if (length(missing_cran)) paste(missing_cran, collapse = ", ") else "none", "\n", sep = "")
cat("Missing Bioconductor packages: ",
    if (length(missing_bioc)) paste(missing_bioc, collapse = ", ") else "none", "\n", sep = "")

# TDCHARM_CHECK_ONLY=1: report only, install nothing.
if (identical(Sys.getenv("TDCHARM_CHECK_ONLY"), "1")) {
    if (length(c(missing_cran, missing_bioc))) quit(status = 1)
    cat("\nAll required packages are installed.\n")
    quit(status = 0)
}

if (length(missing_cran) > 0L) {
    install.packages(missing_cran)
}

if (length(missing_bioc) > 0L && is_missing("BiocManager")) {
    install.packages("BiocManager")
}

if (length(missing_bioc) > 0L) {
    BiocManager::install(
        missing_bioc,
        ask = FALSE,
        update = FALSE
    )
}

still <- c(cran, bioc)[vapply(c(cran, bioc), is_missing, logical(1))]
if (length(still)) {
    stop("Installation failed for: ", paste(still, collapse = ", "))
}

cat("\nPackage check completed.\n")

cat("\nInstalled package versions:\n")
pkgs <- c(cran, bioc)
versions <- vapply(
    pkgs,
    function(p) as.character(utils::packageVersion(p)),
    character(1)
)
print(versions)
