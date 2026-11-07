# Run the revision analyses (scripts 19-27), assuming run_preprocessing.R
# and run_analysis.R have completed in the same working directory.
#
# Usage:
#   GSE303006_WORKDIR=/path/to/analysis Rscript run_revision.R
#
# Optional environment variables (see each script header):
#   TDCHARM_RUN_GSEA=0          skip GO-BP GSEA in 21, 24, 25
#   TDCHARM_RECOMPUTE=1         recompute cached subset runs in 21
#   TDCHARM_EXTERNAL_LINKS=...  external E-P link tables for 23
#   TDCHARM_CHARM_SUPP=...      CHARM Supplementary Tables workbook for 23b
#                               (default: 41586_2026_10322_MOESM3_ESM.xlsx in the
#                               working directory; 23b is skipped if absent)
#   TDCHARM_GROUP_A / _B        exact cell-type labels of the contrast
#
# Script 21 re-runs Stage 1 and Stage 2 (05-08) three times and is by
# far the most expensive step.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
    normalizePath(".")
}
Sys.setenv(TDCHARM_REPO = repo)

workdir <- Sys.getenv("GSE303006_WORKDIR", unset = getwd())
setwd(workdir)

steps <- c(
    "19_stage1_stage2_energy_report.R",
    "20_near_degenerate_rotation.R",
    "21_replicate_subsets.R",
    "21b_factor_matching.R",
    "21c_core12_by_position.R",
    "22_celltype_replicate_joint.R",
    "22b_contrast_by_replicate.R",
    "23_EP_link_benchmark.R",
    "23b_CHARM_link_concordance.R",
    "23c_CHARM_concordance_figure.R",
    "24_modality_ablation.R",
    "24b_ablation_factor_matching.R",
    "25_EP_factor_recurrence.R",
    "26_celltype_component_heatmap.R",
    "27_method_schematic.R"
)

for (s in steps) {
    if (s == "23b_CHARM_link_concordance.R" &&
        !file.exists(Sys.getenv("TDCHARM_CHARM_SUPP", unset = "41586_2026_10322_MOESM3_ESM.xlsx"))) {
        cat("\nSkipping ", s, ": CHARM Supplementary Tables workbook not found ",
            "(download 41586_2026_10322_MOESM3_ESM.xlsx from the article page).\n", sep = "")
        next
    }
    if (s == "23c_CHARM_concordance_figure.R" &&
        !file.exists("GSE303006_M23_rev23b_CHARM_linked_pairs.tsv.gz")) {
        cat("\nSkipping ", s, ": outputs of 23b not found.\n", sep = "")
        next
    }
    cat("\n============================================================\n")
    cat("Running ", s, "\n", sep = "")
    cat("============================================================\n")
    source(file.path(repo, "R", s), local = new.env(parent = globalenv()))
}
