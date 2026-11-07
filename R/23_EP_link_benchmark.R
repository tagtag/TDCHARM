# ============================================================
# 23_EP_link_benchmark.R
# GSE303006 / GENCODE M23
#
# Revision analysis (Reviewer 1, Comment 5).
#
# Are the E-P pairs prioritized by core #12 enriched for supported
# enhancer-gene links, beyond what their genomic distance explains?
#
# Prioritization
#   E-P score of core #12, as in 11:  s_p = CORE[a,b,k,r] U_EP[p,a]
#   Rankings: positive side (s_p), negative side (-s_p), |s_p|.
#
# Reference ("supported") link sets
#   1. ABC-like links computed from this data set (pseudobulk):
#        activity  A_e   = sqrt(mean ATAC CP10K x mean H3K27ac CP10K)
#        contact   C_b   = mean over cells with an observed distance
#                          of distance^(-GAMMA) for the 20-kb bin pair
#        ABC_p           = A_e C_b / sum_{p' with same gene} A_e' C_b'
#      computed for all cells and separately for the two contrasted
#      cell types (Inh Ndnf/Lamp5, Ex L3/4 IT); a pair is supported
#      if ABC_p >= ABC_THRESHOLD. These links come from the same data
#      and are therefore not independent of the decomposition input;
#      they test whether core #12 prioritizes pairs that a standard
#      enhancer-gene model would also call.
#   2. External link tables (optional): any number of tab-separated
#      files given in the environment variable
#         TDCHARM_EXTERNAL_LINKS="file1.tsv,file2.tsv"
#      Each file needs the columns  chr, start, end, gene
#      (enhancer interval in mm10, gene symbol). A candidate pair is
#      supported if its enhancer overlaps a listed interval for the
#      same gene. This is how the E-P links reported in the original
#      CHARM study, or other public mouse-brain link sets, are used.
#
# Statistics, for every reference set and ranking
#   - top N pairs (N = 1000, 5000) versus the rest:
#     Cochran-Mantel-Haenszel test stratified by deciles of
#     |enhancer-TSS distance| (common odds ratio, P)
#   - average precision (AUPRC) of the full ranking, compared with a
#     null obtained by permuting the scores within distance deciles
#     (N_PERM permutations; empirical P), and with the AUPRC of
#     distance alone (closer = higher)
#
# Inputs (working directory): outputs of 03, 05, 07-09, 10
# Outputs:
#   GSE303006_M23_rev23_ABC_like_scores.tsv.gz
#   GSE303006_M23_rev23_reference_sets.tsv
#   GSE303006_M23_rev23_CMH_topN.tsv
#   GSE303006_M23_rev23_AUPRC.tsv
#   GSE303006_M23_rev23_log.txt
# ============================================================

library(data.table)
library(Matrix)

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
ABC_THRESHOLD <- as.numeric(Sys.getenv("TDCHARM_ABC_THRESHOLD", unset = "0.02"))
GAMMA <- 1
SCALE_FACTOR <- 10000
BLOCK_3D <- 2000L
SEED <- 20261103L

EP_FILE <- "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"
CELL_FILE <- "GSE303006_M23_brain_final_4258_cells.rds"
ATAC_FILE <- "GSE303006_M23_brain_ATAC_enhancer_by_cell_AH3D.rds"
H3_FILE <- "GSE303006_M23_brain_H3K27ac_enhancer_by_cell_AH3D.rds"
H5_3D <- "GSE303006_M23_brain_3D_binpair_by_cell.h5"
META_FILE <- "GSE303006_M23_final4258_metadata.tsv.gz"
MAPPING_FILE <- "GSE303006_M23_stage2_EP_mapping_degrees.rds"
CORE_FILE <- "GSE303006_M23_stage2_Tucker_core.rds"
EP_FACTOR_FILE <- "GSE303006_M23_stage2_Tucker_EP_factor.rds"
TOP_CORE_FILE <- "GSE303006_M23_stage2_top100_core.tsv"

OUT_ABC <- "GSE303006_M23_rev23_ABC_like_scores.tsv.gz"
OUT_REF <- "GSE303006_M23_rev23_reference_sets.tsv"
OUT_CMH <- "GSE303006_M23_rev23_CMH_topN.tsv"
OUT_AUPRC <- "GSE303006_M23_rev23_AUPRC.tsv"
OUT_LOG <- "GSE303006_M23_rev23_log.txt"

need <- c(EP_FILE, CELL_FILE, ATAC_FILE, H3_FILE, H5_3D, META_FILE,
          MAPPING_FILE, CORE_FILE, EP_FACTOR_FILE, TOP_CORE_FILE)
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "23_EP_link_benchmark.R")
cat(sprintf("Core rank %d; top N = %s; distance bins = %d; permutations = %d; ABC threshold = %.3f; gamma = %g\n\n",
            CORE_RANK, paste(TOP_N, collapse = ","), N_DIST_BINS, N_PERM, ABC_THRESHOLD, GAMMA))

# ------------------------------------------------------------
# 1. E-P table and core #12 scores
# ------------------------------------------------------------

ep <- as.data.table(readRDS(EP_FILE))
if (!"abs_distance" %in% names(ep)) {
    stopifnot(all(c("enhancer_center", "tss") %in% names(ep)))
    ep[, abs_distance := abs(enhancer_center - tss)]
}
stopifnot(all(c("enhancer_chr", "enhancer_start", "enhancer_end") %in% names(ep)) ||
              !nzchar(Sys.getenv("TDCHARM_EXTERNAL_LINKS")))
map <- as.data.table(readRDS(MAPPING_FILE))
stopifnot(identical(ep$enhancer_id, map$enhancer_id),
          identical(ep$gene_name, map$gene_name))

top <- fread(TOP_CORE_FILE)
cr <- top[core_rank == CORE_RANK]
CORE <- readRDS(CORE_FILE)
UEP <- readRDS(EP_FACTOR_FILE)
G <- CORE[cr$a_EP, cr$b_CELL, cr$k, cr$r_MODAL]
ep[, s := G * UEP[, cr$a_EP]]
cat(sprintf("Core #%d = (EP%02d, CELL%02d, K%02d, MOD%02d), value %.5f\n",
            CORE_RANK, cr$a_EP, cr$b_CELL, cr$k, cr$r_MODAL, G))
rm(UEP); invisible(gc())

ep[, dist_bin := cut(abs_distance,
                     breaks = unique(quantile(abs_distance, seq(0, 1, length.out = N_DIST_BINS + 1))),
                     include.lowest = TRUE, labels = FALSE)]

# ------------------------------------------------------------
# 2. ABC-like scores from this data set
# ------------------------------------------------------------

cells <- as.character(readRDS(CELL_FILE))
meta <- fread(META_FILE)[match(cells, cell)]
groups <- rev_contrast_groups(meta$cell_type)
cell_sets <- list(all = cells,
                  groupA = cells[meta$cell_type %in% groups$A],
                  groupB = cells[meta$cell_type %in% groups$B])
cat("Cells: all = ", length(cell_sets$all), ", A (", paste(groups$A, collapse = ";"), ") = ",
    length(cell_sets$groupA), ", B (", paste(groups$B, collapse = ";"), ") = ",
    length(cell_sets$groupB), "\n", sep = "")

enh <- sort(unique(ep$enhancer_id))
mean_cp10k <- function(file) {
    X <- readRDS(file)
    depth <- Matrix::colSums(X[, cells, drop = FALSE])   # as in 05
    Y <- X[enh, cells, drop = FALSE] %*% Diagonal(x = SCALE_FACTOR / depth)
    colnames(Y) <- cells
    sapply(cell_sets, function(cs) Matrix::rowMeans(Y[, cs, drop = FALSE]))
}
mA <- mean_cp10k(ATAC_FILE)
mH <- mean_cp10k(H3_FILE)
activity <- sqrt(mA * mH)
rownames(activity) <- enh

cat("Reading 3D distances for contact estimates...\n")
bp <- sort(unique(as.integer(ep$binpair_id)))
h5_bp <- as.integer(rhdf5::h5read(H5_3D, "binpair_id"))
h5_cells <- as.character(rhdf5::h5read(H5_3D, "cell"))
rows <- match(bp, h5_bp)
stopifnot(!anyNA(rows))
contact <- matrix(NA_real_, length(bp), length(cell_sets),
                  dimnames = list(as.character(bp), names(cell_sets)))
for (gname in names(cell_sets)) {
    cols <- match(cell_sets[[gname]], h5_cells)
    stopifnot(!anyNA(cols))
    for (ii in split(seq_along(rows), ceiling(seq_along(rows) / BLOCK_3D))) {
        d <- rhdf5::h5read(H5_3D, "distance", index = list(rows[ii], cols))
        d[!is.finite(d) | d <= 0] <- NA_real_
        contact[ii, gname] <- rowMeans(d^(-GAMMA), na.rm = TRUE)
    }
}
rhdf5::h5closeAll()

abc <- ep[, .(enhancer_id, gene_name, binpair_id)]
for (gname in names(cell_sets)) {
    num <- activity[ep$enhancer_id, gname] * contact[as.character(ep$binpair_id), gname]
    num[!is.finite(num)] <- 0
    abc[, (paste0("ABC_", gname)) := num / ave(num, ep$gene_name, FUN = sum)]
}
for (v in grep("^ABC_", names(abc), value = TRUE)) set(abc, which(!is.finite(abc[[v]])), v, 0)
fwrite(abc, OUT_ABC, sep = "\t")

ref_sets <- list(
    ABC_all = abc$ABC_all >= ABC_THRESHOLD,
    ABC_groupA_or_B = abc$ABC_groupA >= ABC_THRESHOLD | abc$ABC_groupB >= ABC_THRESHOLD)

# ------------------------------------------------------------
# 3. External link tables (optional)
# ------------------------------------------------------------

ext <- Sys.getenv("TDCHARM_EXTERNAL_LINKS", unset = "")
if (nzchar(ext)) {
    q <- ep[, .(row = .I, chr = enhancer_chr, start = enhancer_start,
                end = enhancer_end, gene = toupper(gene_name))]
    setkey(q, chr, start, end)
    for (f in trimws(strsplit(ext, ",", fixed = TRUE)[[1]])) {
        L <- fread(f)
        stopifnot(all(c("chr", "start", "end", "gene") %in% names(L)))
        L <- L[, .(chr = as.character(chr), start = as.integer(start),
                   end = as.integer(end), gene = toupper(gene))]
        setkey(L, chr, start, end)
        ov <- foverlaps(q, L, nomatch = 0L)
        hit <- unique(ov[gene == i.gene, row])
        lab <- paste0("EXT_", sub("\\.[^.]*$", "", basename(f)))
        ref_sets[[lab]] <- seq_len(nrow(ep)) %in% hit
        cat("External link set ", lab, ": ", nrow(L), " links, ",
            length(hit), " candidate pairs supported\n", sep = "")
    }
} else {
    cat("No external link table given (TDCHARM_EXTERNAL_LINKS unset).\n")
}

ref_dt <- rbindlist(lapply(names(ref_sets), function(nm) {
    y <- ref_sets[[nm]]
    data.table(reference = nm, n_supported = sum(y), fraction = mean(y),
               median_distance_supported = median(ep$abs_distance[y]),
               median_distance_other = median(ep$abs_distance[!y]))
}))
cat("\nReference sets:\n")
print(ref_dt)

# ------------------------------------------------------------
# 4. Statistics
# ------------------------------------------------------------

avg_precision <- function(score, y) {
    o <- order(score, decreasing = TRUE)
    yy <- y[o]
    tp <- cumsum(yy)
    sum((tp / seq_along(yy))[yy]) / sum(yy)
}

rankings <- list(positive = ep$s, negative = -ep$s, absolute = abs(ep$s))
bins <- ep$dist_bin
by_bin <- split(seq_len(nrow(ep)), bins)

cmh_rows <- list(); ap_rows <- list()
set.seed(SEED)
for (rn in names(ref_sets)) {
    y <- ref_sets[[rn]]
    if (sum(y) < 10L || sum(!y) < 10L) {
        cat("Skipping ", rn, ": fewer than 10 supported or unsupported pairs\n", sep = "")
        next
    }
    ap_dist <- avg_precision(-ep$abs_distance, y)
    for (sd in names(rankings)) {
        sc <- rankings[[sd]]
        for (n in TOP_N) {
            if (n >= length(sc)) next
            topf <- seq_along(sc) %in% order(sc, decreasing = TRUE)[seq_len(n)]
            arr <- table(factor(topf, c(TRUE, FALSE)), factor(y, c(TRUE, FALSE)),
                         factor(bins))
            arr <- arr[, , apply(arr, 3, function(z) all(rowSums(z) > 0) && all(colSums(z) > 0)), drop = FALSE]
            mh <- mantelhaen.test(arr, exact = FALSE)
            cmh_rows[[length(cmh_rows) + 1L]] <- data.table(
                reference = rn, ranking = sd, top_n = n,
                supported_in_top = sum(y[topf]), fraction_in_top = mean(y[topf]),
                fraction_overall = mean(y),
                common_OR = as.numeric(mh$estimate),
                OR_CI_low = mh$conf.int[1], OR_CI_high = mh$conf.int[2],
                CMH_p = mh$p.value)
        }
        ap <- avg_precision(sc, y)
        null <- numeric(N_PERM)
        for (i in seq_len(N_PERM)) {
            sp <- sc
            for (ii in by_bin) sp[ii] <- sc[ii][sample.int(length(ii))]
            null[i] <- avg_precision(sp, y)
        }
        ap_rows[[length(ap_rows) + 1L]] <- data.table(
            reference = rn, ranking = sd, AUPRC = ap,
            AUPRC_null_distance_matched_mean = mean(null),
            AUPRC_null_sd = sd(null),
            perm_p = (1 + sum(null >= ap)) / (1 + N_PERM),
            AUPRC_distance_only = ap_dist, baseline_prevalence = mean(y))
    }
}
cmh <- rbindlist(cmh_rows)
aup <- rbindlist(ap_rows)
if (!nrow(cmh)) cmh <- data.table(reference = character(), ranking = character(), top_n = integer(),
                                  supported_in_top = integer(), fraction_in_top = numeric(),
                                  fraction_overall = numeric(), common_OR = numeric(), CMH_p = numeric())

cat("\nTop-N enrichment, stratified by distance decile (CMH):\n")
print(cmh[, .(reference, ranking, top_n, supported_in_top,
              fraction_in_top = round(fraction_in_top, 4),
              fraction_overall = round(fraction_overall, 4),
              OR = round(common_OR, 3), CMH_p = signif(CMH_p, 3))])
cat("\nAUPRC versus distance-matched permutation null:\n")
print(aup[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)])

fwrite(ref_dt, OUT_REF, sep = "\t")
fwrite(cmh, OUT_CMH, sep = "\t")
fwrite(aup, OUT_AUPRC, sep = "\t")
cat("\nSaved outputs with prefix GSE303006_M23_rev23_.\n")
rev_end_log(log_con)
