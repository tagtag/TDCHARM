# Run the hierarchical tensor analysis and manuscript-level
# downstream analyses, assuming preprocessing stages 01--04 are complete.
#
# Usage:
#   GSE303006_WORKDIR=/path/to/analysis Rscript run_analysis.R

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
    "05_stage1_independent_K20.R",
    "06_sign_only_stage1_K20.R",
    "07_prepare_stage2_implicit.R",
    "08_stage2_implicit_HOSVD.R",
    "09_trace_core_to_original.R",
    "10_celltype_replicate_annotation.R",
    "11_core12_allEP_gene_rank.R",
    "12_core12_mouse_GO_BP_GSEA.R",
    "13_core12_GO_reduce_redundancy.R",
    "14_core12_GO_process_gene_overlap.R",
    "15_core12_original4modality_backprojection.R",
    "16_core12_replicate_reproducibility.R",
    "17_top100_CELLcomponent_effects.R",
    "18_top100_publication_figure.R"
)

for (s in steps) {
    cat("\n============================================================\n")
    cat("Running ", s, "\n", sep = "")
    cat("============================================================\n")
    source(file.path(repo, "R", s), chdir = FALSE)
}
