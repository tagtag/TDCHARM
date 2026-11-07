# ============================================================
# 23b_CHARM_link_concordance.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 1, Comment 5), complement to 23.
#
# Concordance between the E-P pairs prioritized by core #12 and the
# enhancer-gene links reported in the original CHARM study
# (Chen et al., Nature 2026, doi:10.1038/s41586-026-10322-z,
# Supplementary Table 7: 4,954 unique predicted E-P pairs, 5-kb bins,
# with presence/absence in 16 cell types).
#
# The CHARM links are themselves predictions from the same data set
# (Poisson regression + SHAP, threshold calibrated on 180 mESC pairs,
# Supplementary Table 6), so this is a concordance analysis between
# two analyses of GSE303006, not a validation against ground truth.
#
# The supplementary workbook is not redistributed with this
# repository. Download it from the article page
# ("Supplementary Tables", file 41586_2026_10322_MOESM3_ESM.xlsx)
# into the working directory, or give its path in
#    TDCHARM_CHARM_SUPP=/path/to/41586_2026_10322_MOESM3_ESM.xlsx
#
# Matching
#   A candidate pair (enhancer e, gene g) of this study is CHARM-linked
#   if e overlaps a CHARM 5-kb bin [start+1, start+5000] listed for the
#   same gene symbol (case-insensitive).
#
# Universes (CHARM tested only its 2,000 most variable genes within
# +-2 Mb; this study uses all genes within 1 Mb)
#   all_pairs     all 730,969 candidate pairs
#   charm_genes   pairs whose gene has at least one link in Table 7
#
# Reference sets
#   CHARM_any            linked in any cell type
#   CHARM_intergenic     linked in any cell type, bin outside gene body
#   CHARM_groupA         linked in the group-A cell type(s) (Ndnf/Lamp5)
#   CHARM_groupB         linked in the group-B cell type(s) (L3/4 IT)
#   CHARM_groupA_or_B
#
# Statistics (as in 23), for each universe x reference set x ranking
# (positive s_p, negative -s_p, |s_p|):
#   - top N pairs vs rest, CMH test stratified by distance decile
#   - AUPRC vs distance-matched permutation null and vs distance only
# and, among CHARM-linked pairs, core #12 score of links specific to
# group A versus links specific to group B (Wilcoxon, rank-biserial).
#
# Inputs (working directory): outputs of 02, 07-09 and the workbook
# Outputs:
#   GSE303006_M23_rev23b_CHARM_links_parsed.tsv.gz
#   GSE303006_M23_rev23b_CHARM_coverage.tsv
#   GSE303006_M23_rev23b_reference_sets.tsv
#   GSE303006_M23_rev23b_CMH_topN.tsv
#   GSE303006_M23_rev23b_AUPRC.tsv
#   GSE303006_M23_rev23b_groupA_vs_groupB_links.tsv
#   GSE303006_M23_rev23b_CHARM_linked_pairs.tsv.gz (input of 23c)
#   GSE303006_M23_rev23b_log.txt
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

# ------------------------------------------------------------
# 0. Settings / files
# ------------------------------------------------------------

CORE_RANK <- 12L
TOP_N <- as.integer(strsplit(Sys.getenv("TDCHARM_TOP_N", unset = "1000,5000"), ",")[[1]])
N_DIST_BINS <- 10L
N_PERM <- as.integer(Sys.getenv("TDCHARM_N_PERM", unset = "200"))
BIN_SIZE <- 5000L
SEED <- 20261104L

SUPP_FILE <- Sys.getenv("TDCHARM_CHARM_SUPP", unset = "41586_2026_10322_MOESM3_ESM.xlsx")
SUPP_SHEET <- "Supplementary Table 7"

EP_FILE <- "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
EP_FACTOR_FILE <- "GSE303006_M23_stage2_Tucker_EP_factor.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"

OUT_PARSED <- "GSE303006_M23_rev23b_CHARM_links_parsed.tsv.gz"
OUT_COV <- "GSE303006_M23_rev23b_CHARM_coverage.tsv"
OUT_REF <- "GSE303006_M23_rev23b_reference_sets.tsv"
OUT_CMH <- "GSE303006_M23_rev23b_CMH_topN.tsv"
OUT_AUPRC <- "GSE303006_M23_rev23b_AUPRC.tsv"
OUT_AB <- "GSE303006_M23_rev23b_groupA_vs_groupB_links.tsv"
OUT_PAIRS <- "GSE303006_M23_rev23b_CHARM_linked_pairs.tsv.gz"
OUT_LOG <- "GSE303006_M23_rev23b_log.txt"

need <- c(SUPP_FILE, EP_FILE, CORE_FILE, EP_FACTOR_FILE, TOP_CORE_FILE, META_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "),
                                  "\n(the CHARM workbook must be downloaded from the article page; see header)")
if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required (see 00_install_packages.R)")

log_con <- rev_start_log(OUT_LOG, "23b_CHARM_link_concordance.R")
cat(sprintf("Core rank %d; top N = %s; distance bins = %d; permutations = %d\n",
            CORE_RANK, paste(TOP_N, collapse = ","), N_DIST_BINS, N_PERM))
cat("CHARM workbook: ", normalizePath(SUPP_FILE), "\n", sep = "")
cat("CHARM workbook md5: ", unname(tools::md5sum(SUPP_FILE)), "\n\n", sep = "")

# ------------------------------------------------------------
# 1. Parse Supplementary Table 7
#    row 1 title; row 2: Chrom, Start, End, Within gene-body, Gene,
#    "Predicted to be linked in specific cell type"; row 3: cell types
# ------------------------------------------------------------

hdr <- as.data.frame(readxl::read_excel(SUPP_FILE, sheet = SUPP_SHEET, col_names = FALSE,
                                        n_max = 3, .name_repair = "minimal"))
raw <- as.data.frame(readxl::read_excel(SUPP_FILE, sheet = SUPP_SHEET, col_names = FALSE,
                                        skip = 3, .name_repair = "minimal"))
stopifnot(identical(as.character(unlist(hdr[2, 1:5])),
                    c("Chrom", "Start", "End", "Within gene-body", "Gene")))
ct_cols <- 6:ncol(hdr)
cell_types <- as.character(unlist(hdr[3, ct_cols]))
stopifnot(!anyNA(cell_types))

charm <- data.table(chr = as.character(raw[[1]]),
                    bin_start = as.integer(round(as.numeric(raw[[2]]))),
                    within_gene_body = as.logical(raw[[4]]),
                    gene = as.character(raw[[5]]))
# End is stored as a formula (=Start+5000); recompute it
charm[, bin_end := bin_start + BIN_SIZE]
flags <- as.data.table(lapply(raw[ct_cols], function(x) as.integer(as.numeric(x))))
setnames(flags, cell_types)
charm <- cbind(charm, flags)
charm <- charm[!is.na(chr) & !is.na(bin_start) & !is.na(gene)]
charm[, n_cell_types := rowSums(.SD), .SDcols = cell_types]
stopifnot(all(charm$n_cell_types >= 1L))
cat(sprintf("CHARM Table 7: %d unique E-P pairs, %d genes, %d cell-type linkages, %d cell types\n",
            nrow(charm), uniqueN(charm$gene), sum(charm$n_cell_types), length(cell_types)))
print(sort(colSums(charm[, ..cell_types]), decreasing = TRUE))
fwrite(charm, OUT_PARSED, sep = "\t")

# Contrast groups, same rule as the rest of the revision scripts
meta <- fread(META_FILE)
groups <- rev_contrast_groups(meta$cell_type)
# Labels are compared after removing non-alphanumeric characters, and a
# CHARM column matches if this study's label ends with it
# (e.g. "Ex_L3/4_IT" -> "L3/4 IT", "Inh_Ndnf/Lamp5" -> "Ndnf/Lamp5").
norm_lab <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
match_ct <- function(g) {
    m <- cell_types[vapply(cell_types, function(ct)
        any(endsWith(norm_lab(g), norm_lab(ct))), logical(1))]
    if (!length(m)) stop("No CHARM cell-type column matches: ", paste(g, collapse = "; "))
    m
}
ctA <- match_ct(groups$A); ctB <- match_ct(groups$B)
cat("\nGroup A: ", paste(groups$A, collapse = ";"), " -> CHARM column(s) ", paste(ctA, collapse = ";"), "\n",
    "Group B: ", paste(groups$B, collapse = ";"), " -> CHARM column(s) ", paste(ctB, collapse = ";"), "\n", sep = "")
charm[, inA := rowSums(.SD) > 0, .SDcols = ctA]
charm[, inB := rowSums(.SD) > 0, .SDcols = ctB]

# ------------------------------------------------------------
# 2. Candidate pairs of this study and core #12 scores
# ------------------------------------------------------------

ep <- as.data.table(readRDS(EP_FILE))
stopifnot(all(c("enhancer_chr", "enhancer_start", "enhancer_end", "gene_name") %in% names(ep)))
if (!"abs_distance" %in% names(ep)) ep[, abs_distance := abs(enhancer_center - tss)]
top <- fread(TOP_CORE_FILE)
cr <- top[core_rank == CORE_RANK]
CORE <- readRDS(CORE_FILE)
UEP <- readRDS(EP_FACTOR_FILE)
G <- CORE[cr$a_EP, cr$b_CELL, cr$k, cr$r_MODAL]
ep[, s := G * UEP[, cr$a_EP]]
rm(UEP); invisible(gc())
cat(sprintf("\nCore #%d = (EP%02d, CELL%02d, K%02d, MOD%02d), value %.5f; %d candidate pairs, %d genes\n",
            CORE_RANK, cr$a_EP, cr$b_CELL, cr$k, cr$r_MODAL, G, nrow(ep), uniqueN(ep$gene_name)))

# ------------------------------------------------------------
# 3. Overlap
# ------------------------------------------------------------

q <- ep[, .(row = .I, chr = as.character(enhancer_chr), start = as.integer(enhancer_start),
            end = as.integer(enhancer_end), G = toupper(gene_name))]
L <- charm[, .(cid = .I, chr, start = bin_start + 1L, end = bin_end, G = toupper(gene))]
setkey(L, chr, start, end)
ov <- foverlaps(q, L, nomatch = 0L)[G == i.G, .(row, cid)]

ep[, `:=`(CHARM_any = FALSE, CHARM_intergenic = FALSE, CHARM_groupA = FALSE, CHARM_groupB = FALSE)]
if (nrow(ov)) {
    agg <- ov[, .(any = TRUE,
                  intergenic = any(!charm$within_gene_body[cid]),
                  A = any(charm$inA[cid]), B = any(charm$inB[cid])), by = row]
    ep[agg$row, `:=`(CHARM_any = agg$any, CHARM_intergenic = agg$intergenic,
                     CHARM_groupA = agg$A, CHARM_groupB = agg$B)]
}
ep[, CHARM_groupA_or_B := CHARM_groupA | CHARM_groupB]

# Coverage of the CHARM links by the candidate set of this study
gene_tss <- unique(ep[, .(G = toupper(gene_name), chr = as.character(enhancer_chr), tss)])
gene_tss <- gene_tss[!duplicated(G)]
charm[, G := toupper(gene)]
charm[, gene_in_study := G %in% gene_tss$G]
charm[gene_tss, tss := i.tss, on = "G"]
charm[, bin_tss_distance := abs(bin_start + BIN_SIZE / 2 - tss)]
charm[, recovered := seq_len(.N) %in% unique(ov$cid)]
cov <- rbind(
    charm[, .(subset = "all", n = .N, gene_in_study = sum(gene_in_study),
              within_1Mb = sum(bin_tss_distance <= 1e6, na.rm = TRUE), recovered = sum(recovered))],
    charm[inA == TRUE, .(subset = "groupA", n = .N, gene_in_study = sum(gene_in_study),
                         within_1Mb = sum(bin_tss_distance <= 1e6, na.rm = TRUE), recovered = sum(recovered))],
    charm[inB == TRUE, .(subset = "groupB", n = .N, gene_in_study = sum(gene_in_study),
                         within_1Mb = sum(bin_tss_distance <= 1e6, na.rm = TRUE), recovered = sum(recovered))])
cat("\nCoverage of CHARM Table 7 pairs by this study's candidate E-P pairs\n",
    "(gene_in_study: gene among this study's genes; within_1Mb: bin centre within 1 Mb of\n",
    " this study's TSS; recovered: bin overlaps a candidate enhancer of the same gene)\n", sep = "")
print(cov)
fwrite(cov, OUT_COV, sep = "\t")

# ------------------------------------------------------------
# 4. Enrichment statistics (same procedures as 23)
# ------------------------------------------------------------

avg_precision <- function(score, y) {
    o <- order(score, decreasing = TRUE)
    yy <- y[o]
    tp <- cumsum(yy)
    sum((tp / seq_along(yy))[yy]) / sum(yy)
}

ref_names <- c("CHARM_any", "CHARM_intergenic", "CHARM_groupA", "CHARM_groupB", "CHARM_groupA_or_B")
charm_genes <- unique(charm$G)
universes <- list(all_pairs = seq_len(nrow(ep)),
                  charm_genes = which(toupper(ep$gene_name) %in% charm_genes))

ref_rows <- list(); cmh_rows <- list(); ap_rows <- list()
set.seed(SEED)
for (un in names(universes)) {
    idx <- universes[[un]]
    E <- ep[idx]
    E[, dist_bin := cut(abs_distance,
                        breaks = unique(quantile(abs_distance, seq(0, 1, length.out = N_DIST_BINS + 1))),
                        include.lowest = TRUE, labels = FALSE)]
    bins <- E$dist_bin
    by_bin <- split(seq_len(nrow(E)), bins)
    rankings <- list(positive = E$s, negative = -E$s, absolute = abs(E$s))
    cat(sprintf("\nUniverse %s: %d pairs, %d genes\n", un, nrow(E), uniqueN(E$gene_name)))
    for (rn in ref_names) {
        y <- E[[rn]]
        ref_rows[[length(ref_rows) + 1L]] <- data.table(
            universe = un, reference = rn, n_pairs = nrow(E), n_supported = sum(y), fraction = mean(y),
            median_distance_supported = if (any(y)) as.numeric(median(E$abs_distance[y])) else NA_real_,
            median_distance_other = as.numeric(median(E$abs_distance[!y])))
        if (sum(y) < 10L || sum(!y) < 10L) {
            cat("  Skipping ", rn, ": fewer than 10 supported or unsupported pairs\n", sep = "")
            next
        }
        ap_dist <- avg_precision(-E$abs_distance, y)
        for (sd in names(rankings)) {
            sc <- rankings[[sd]]
            for (n in TOP_N) {
                if (n >= length(sc)) next
                topf <- seq_along(sc) %in% order(sc, decreasing = TRUE)[seq_len(n)]
                arr <- table(factor(topf, c(TRUE, FALSE)), factor(y, c(TRUE, FALSE)), factor(bins))
                arr <- arr[, , apply(arr, 3, function(z) all(rowSums(z) > 0) && all(colSums(z) > 0)), drop = FALSE]
                mh <- if (dim(arr)[3] >= 1L) mantelhaen.test(arr, exact = FALSE) else NULL
                cmh_rows[[length(cmh_rows) + 1L]] <- data.table(
                    universe = un, reference = rn, ranking = sd, top_n = n,
                    supported_in_top = sum(y[topf]), fraction_in_top = mean(y[topf]),
                    fraction_overall = mean(y),
                    common_OR = if (is.null(mh)) NA_real_ else as.numeric(mh$estimate),
                    OR_CI_low = if (is.null(mh)) NA_real_ else mh$conf.int[1],
                    OR_CI_high = if (is.null(mh)) NA_real_ else mh$conf.int[2],
                    CMH_p = if (is.null(mh)) NA_real_ else mh$p.value)
            }
            ap <- avg_precision(sc, y)
            null <- numeric(N_PERM)
            for (i in seq_len(N_PERM)) {
                sp <- sc
                for (ii in by_bin) sp[ii] <- sc[ii][sample.int(length(ii))]
                null[i] <- avg_precision(sp, y)
            }
            ap_rows[[length(ap_rows) + 1L]] <- data.table(
                universe = un, reference = rn, ranking = sd, AUPRC = ap,
                AUPRC_null_distance_matched_mean = mean(null), AUPRC_null_sd = sd(null),
                perm_p = (1 + sum(null >= ap)) / (1 + N_PERM),
                AUPRC_distance_only = ap_dist, baseline_prevalence = mean(y))
        }
    }
}
ref_dt <- rbindlist(ref_rows)
cmh <- rbindlist(cmh_rows, fill = TRUE)
if (!nrow(cmh)) cmh <- data.table(universe = character(), reference = character(), ranking = character(),
                                  top_n = integer(), supported_in_top = integer(), fraction_in_top = numeric(),
                                  fraction_overall = numeric(), common_OR = numeric(), OR_CI_low = numeric(),
                                  OR_CI_high = numeric(), CMH_p = numeric())
aup <- rbindlist(ap_rows, fill = TRUE)

cat("\nReference sets:\n"); print(ref_dt)
if (nrow(cmh)) {
    cat("\nTop-N enrichment, stratified by distance decile (CMH):\n")
    print(cmh[, .(universe, reference, ranking, top_n, supported_in_top,
                  fraction_in_top = round(fraction_in_top, 4), fraction_overall = round(fraction_overall, 4),
                  OR = round(common_OR, 3), CMH_p = signif(CMH_p, 3))], nrows = 500)
}
if (nrow(aup)) {
    cat("\nAUPRC versus distance-matched permutation null:\n")
    print(aup[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], nrows = 500)
}

# ------------------------------------------------------------
# 5. Among CHARM-linked candidate pairs: group-A-specific vs
#    group-B-specific links, core #12 score
# ------------------------------------------------------------

# percentile of s among all candidate pairs (50 = middle of the ranking)
ep[, s_percentile := 100 * (frank(s, ties.method = "average") - 0.5) / .N]
ab <- ep[CHARM_any == TRUE, .(enhancer_id, gene_name, abs_distance, s, s_percentile,
                              CHARM_intergenic, CHARM_groupA, CHARM_groupB,
                              cls = fifelse(CHARM_groupA & !CHARM_groupB, "A_only",
                                    fifelse(CHARM_groupB & !CHARM_groupA, "B_only",
                                    fifelse(CHARM_groupA & CHARM_groupB, "A_and_B", "other"))))]
ab_sum <- ab[, .(n = .N, mean_s = mean(s), median_s = median(s),
                 frac_positive = mean(s > 0), median_distance = as.numeric(median(abs_distance))), by = cls][order(cls)]
cat("\nCore #12 score of CHARM-linked candidate pairs by CHARM cell-type class:\n"); print(ab_sum)
xa <- ab[cls == "A_only", s]; xb <- ab[cls == "B_only", s]
if (length(xa) >= 3L && length(xb) >= 3L) {
    wt <- wilcox.test(xa, xb)
    ab_test <- data.table(n_A_only = length(xa), n_B_only = length(xb),
                          median_s_A_only = median(xa), median_s_B_only = median(xb),
                          rank_biserial_A_minus_B = rev_rank_biserial(xa, xb),
                          wilcox_p = wt$p.value)
    cat("\nA-only vs B-only links (Wilcoxon rank-sum; rank-biserial > 0: A-only higher):\n"); print(ab_test)
} else {
    ab_test <- data.table(n_A_only = length(xa), n_B_only = length(xb))
    cat("\nToo few A-only / B-only links among candidate pairs for a test.\n")
}
fwrite(ab_test, OUT_AB, sep = "\t")
fwrite(ab, OUT_PAIRS, sep = "\t")
fwrite(ab_sum, sub("\\.tsv$", "_summary.tsv", OUT_AB), sep = "\t")

fwrite(ref_dt, OUT_REF, sep = "\t")
fwrite(cmh, OUT_CMH, sep = "\t")
fwrite(aup, OUT_AUPRC, sep = "\t")
cat("\nSaved outputs with prefix GSE303006_M23_rev23b_.\n")
rev_end_log(log_con)
