# ============================================================
# 22b_contrast_by_replicate.R
# GSE303006 / GENCODE M23
#
# Follow-up to 22 (Reviewer 1, Comment 4).
#
# Script 22 found a replicate effect on the CELL03 loading WITHIN one
# of the two cell types that define the core #12 contrast
# (Inh Ndnf/Lamp5). This script asks how large that within-type
# replicate effect is compared with the between-type contrast itself:
#
#   1. CELL03 loading (oriented: group A higher, as in 21) summarized
#      per cell type of the contrast and per replicate
#      (n, mean, median, IQR)
#   2. within each replicate: rank-biserial effect of group A versus
#      group B and Wilcoxon P (is the contrast present in every
#      replicate, and how strong?)
#   3. within each group: pairwise replicate comparisons
#      (rank-biserial, Wilcoxon P, BH), i.e. which replicate differs
#   4. the replicate shift relative to the between-type difference:
#      max |difference of replicate medians| within a group divided by
#      |median(A) - median(B)| over all cells
#
# The cell component and the cell types are those of manuscript
# core #12 (top-100 table, core rank 12) and the pre-specified
# contrast of 21 (revision_lib.R).
#
# Inputs (working directory): outputs of 08-10.
# Outputs:
#   GSE303006_M23_rev22b_summary_by_type_replicate.tsv
#   GSE303006_M23_rev22b_contrast_within_replicate.tsv
#   GSE303006_M23_rev22b_replicate_pairs_within_type.tsv
#   GSE303006_M23_rev22b_log.txt
# ============================================================

library(data.table)

local({
    lib <- file.path(Sys.getenv("TDCHARM_REPO", unset = "."), "R", "revision_lib.R")
    if (!exists("rev_find_repo")) {
        if (!file.exists(lib)) {
            fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
            if (length(fa)) lib <- file.path(dirname(normalizePath(sub("^--file=", "", fa[1]))), "revision_lib.R")
            if (!file.exists(lib) && exists("repo", envir = globalenv()))
                lib <- file.path(get("repo", envir = globalenv()), "R", "revision_lib.R")
        }
        sys.source(lib, envir = globalenv())
    }
})

CORE_RANK <- 12L
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_SUM <- "GSE303006_M23_rev22b_summary_by_type_replicate.tsv"
OUT_CON <- "GSE303006_M23_rev22b_contrast_within_replicate.tsv"
OUT_PAIR <- "GSE303006_M23_rev22b_replicate_pairs_within_type.tsv"
OUT_LOG <- "GSE303006_M23_rev22b_log.txt"

need <- c(SMALL_FILE, TOP_CORE_FILE, META_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "22b_contrast_by_replicate.R")

SF <- readRDS(SMALL_FILE)
b <- fread(TOP_CORE_FILE)[core_rank == CORE_RANK, b_CELL]
comp <- sprintf("CELL%02d", b)
cells <- rownames(SF$U_CELL)
meta <- fread(META_FILE)[match(cells, cell)]
groups <- rev_contrast_groups(meta$cell_type)
grp <- fifelse(meta$cell_type %in% groups$A, "A",
               fifelse(meta$cell_type %in% groups$B, "B", NA_character_))
x <- SF$U_CELL[, b]
if (mean(x[grp == "A"], na.rm = TRUE) < mean(x[grp == "B"], na.rm = TRUE)) x <- -x

cat("Component: ", comp, " (core rank ", CORE_RANK, "), oriented so that group A is higher\n", sep = "")
cat("Group A: ", paste(groups$A, collapse = "; "), "\nGroup B: ", paste(groups$B, collapse = "; "), "\n\n", sep = "")

d <- data.table(cell = cells, group = grp, replicate = meta$replicate, x = x)[!is.na(group)]
reps <- sort(unique(d$replicate))

# 1. summary
summ <- d[, .(n = .N, mean = mean(x), median = median(x),
              q25 = quantile(x, 0.25), q75 = quantile(x, 0.75)),
          by = .(group, replicate)][order(group, replicate)]
summ <- rbind(summ, d[, .(replicate = "all", n = .N, mean = mean(x), median = median(x),
                          q25 = quantile(x, 0.25), q75 = quantile(x, 0.75)), by = group])
cat("CELL loading by group and replicate:\n")
print(summ[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)])

# 2. contrast within each replicate
con <- rbindlist(lapply(c(reps, "all"), function(r) {
    z <- if (r == "all") d else d[replicate == r]
    xa <- z[group == "A", x]; xb <- z[group == "B", x]
    data.table(replicate = r, nA = length(xa), nB = length(xb),
               rank_biserial_A_vs_B = rev_rank_biserial(xa, xb),
               wilcox_p = suppressWarnings(wilcox.test(xa, xb)$p.value))
}))
cat("\nGroup A versus group B within each replicate:\n")
print(con[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)])

# 3. replicate pairs within each group
pairs <- combn(reps, 2, simplify = FALSE)
pr <- rbindlist(lapply(c("A", "B"), function(g) rbindlist(lapply(pairs, function(p) {
    x1 <- d[group == g & replicate == p[1], x]; x2 <- d[group == g & replicate == p[2], x]
    data.table(group = g, replicate_1 = p[1], replicate_2 = p[2],
               n1 = length(x1), n2 = length(x2),
               median_diff = median(x1) - median(x2),
               rank_biserial = rev_rank_biserial(x1, x2),
               wilcox_p = suppressWarnings(wilcox.test(x1, x2)$p.value))
}))))
pr[, p_BH := p.adjust(wilcox_p, "BH")]
cat("\nPairwise replicate comparisons within each group:\n")
print(pr[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)])

# 4. replicate shift relative to the contrast
between <- abs(summ[group == "A" & replicate == "all", median] -
               summ[group == "B" & replicate == "all", median])
shift <- pr[, .(max_abs_median_diff = max(abs(median_diff))), by = group]
shift[, ratio_to_between_type_difference := max_abs_median_diff / between]
cat(sprintf("\n|median(A) - median(B)| over all cells = %.4g\n", between))
cat("Largest replicate shift within each group, relative to it:\n")
print(shift[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 4) else v)])

fwrite(summ, OUT_SUM, sep = "\t")
fwrite(con, OUT_CON, sep = "\t")
fwrite(merge(pr, shift, by = "group"), OUT_PAIR, sep = "\t")
cat("\nSaved outputs with prefix GSE303006_M23_rev22b_.\n")
rev_end_log(log_con)
