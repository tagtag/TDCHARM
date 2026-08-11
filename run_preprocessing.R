# Run preprocessing stages 01--04.
#
# Usage:
#   GSE303006_WORKDIR=/path/to/analysis Rscript run_preprocessing.R
#
# The working directory should contain the downloaded GSE303006 files,
# GENCODE M23 annotation, and reconstructed 3D directory.
# Peak calling is NOT rerun by default in script 01.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
    normalizePath(".")
}

workdir <- Sys.getenv("GSE303006_WORKDIR", unset = getwd())
setwd(workdir)

steps <- c(
    "01_prepare_annotation_independent_inputs.R",
    "02_build_M23_EP_and_3D_summary.R",
    "03_finalize_EP_RNA_and_cell_QC.R",
    "04_build_3D_HDF5.R"
)

for (s in steps) {
    cat("\n============================================================\n")
    cat("Running ", s, "\n", sep = "")
    cat("============================================================\n")
    source(file.path(repo, "R", s), chdir = FALSE)
}
