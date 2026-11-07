# ============================================================
# 21b_factor_matching.R
# GSE303006 / GENCODE M23
#
# Follow-up to 21_replicate_subsets.R (Reviewer 1, Comment 3).
#
# The index of a Stage-2 factor (EP01, EP02, ...; CELL01, ...) is only
# its rank within one decomposition. A factor of a subset run may be the
# same vector as a full-data factor under a different index. This script
# compares the factor VECTORS directly, without using any selection rule:
#
#   E-P factors  : |Pearson correlation| between every full-data E-P
#                  factor and every subset E-P factor over the 730,969
#                  E-P pairs (same rows in all runs)
#   cell factors : |Pearson correlation| over the cells shared by the
#                  full data and the subset
#
# For every full-data factor the best-matching subset factor and its
# |r| are reported, together with the second-best |r| (a clear match
# has a large best |r| and a small second-best |r|). The factors used by
# the core selected in each run (21) are listed with their matches.
#
# Inputs (working directory): outputs of 08 (full data), the sub-
# directories rev21_* and GSE303006_M23_rev21_selection.tsv from 21.
# Outputs:
#   GSE303006_M23_rev21b_EP_factor_abs_correlation.tsv
#   GSE303006_M23_rev21b_CELL_factor_abs_correlation.tsv
#   GSE303006_M23_rev21b_best_matches.tsv
#   GSE303006_M23_rev21b_log.txt
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

SUBSETS <- c("R1R2", "R1R3", "R2R3")
EP_FACTOR_FILE <- "GSE303006_M23_stage2_Tucker_EP_factor.rds"
SMALL_FILE <- "GSE303006_M23_stage2_Tucker_small_factors.rds"
SEL_FILE <- "GSE303006_M23_rev21_selection.tsv"

OUT_EP <- "GSE303006_M23_rev21b_EP_factor_abs_correlation.tsv"
OUT_CELL <- "GSE303006_M23_rev21b_CELL_factor_abs_correlation.tsv"
OUT_BEST <- "GSE303006_M23_rev21b_best_matches.tsv"
OUT_LOG <- "GSE303006_M23_rev21b_log.txt"

need <- c(EP_FACTOR_FILE, SMALL_FILE,
          file.path(paste0("rev21_", SUBSETS), EP_FACTOR_FILE),
          file.path(paste0("rev21_", SUBSETS), SMALL_FILE))
if (any(!file.exists(need))) stop("Missing: ", paste(need[!file.exists(need)], collapse = ", "))

log_con <- rev_start_log(OUT_LOG, "21b_factor_matching.R")

abs_cor <- function(A, B) abs(cor(A, B))

best_table <- function(R, mode, subset) {
    rbindlist(lapply(seq_len(nrow(R)), function(i) {
        o <- order(R[i, ], decreasing = TRUE)
        data.table(mode = mode, subset = subset,
                   full_factor = rownames(R)[i],
                   best_subset_factor = colnames(R)[o[1]],
                   best_abs_r = R[i, o[1]],
                   second_subset_factor = colnames(R)[o[2]],
                   second_abs_r = R[i, o[2]])
    }))
}

UEP_full <- readRDS(EP_FACTOR_FILE)
SF_full <- readRDS(SMALL_FILE)
UC_full <- SF_full$U_CELL

ep_rows <- list(); cell_rows <- list(); best_rows <- list()
for (s in SUBSETS) {
    d <- paste0("rev21_", s)
    UEP_s <- readRDS(file.path(d, EP_FACTOR_FILE))
    stopifnot(nrow(UEP_s) == nrow(UEP_full))
    R_ep <- abs_cor(UEP_full, UEP_s)
    dimnames(R_ep) <- list(colnames(UEP_full), colnames(UEP_s))

    UC_s <- readRDS(file.path(d, SMALL_FILE))$U_CELL
    common <- intersect(rownames(UC_full), rownames(UC_s))
    R_cell <- abs_cor(UC_full[common, , drop = FALSE], UC_s[common, , drop = FALSE])
    dimnames(R_cell) <- list(colnames(UC_full), colnames(UC_s))

    ep_rows[[s]] <- cbind(subset = s, as.data.table(R_ep, keep.rownames = "full_factor"))
    cell_rows[[s]] <- cbind(subset = s, as.data.table(R_cell, keep.rownames = "full_factor"))
    best_rows[[s]] <- rbind(best_table(R_ep, "EP", s), best_table(R_cell, "CELL", s))
    rm(UEP_s); invisible(gc())
}
best <- rbindlist(best_rows)

cat("Best-matching subset factor for every full-data factor (|r|):\n\n")
for (m in c("EP", "CELL")) {
    w <- dcast(best[mode == m], full_factor ~ subset,
               value.var = c("best_subset_factor", "best_abs_r"))
    num <- grep("best_abs_r", names(w), value = TRUE)
    w[, (num) := lapply(.SD, round, 3), .SDcols = num]
    cat(m, "factors:\n"); print(w[]); cat("\n")
}

if (file.exists(SEL_FILE)) {
    sel <- fread(SEL_FILE)
    cat("Factors used by the core selected in each run (21), and their matches:\n")
    f <- sel[run == "FULL"]
    cat(sprintf("  FULL: EP%02d, CELL%02d\n", f$EP, f$CELL))
    for (s in SUBSETS) {
        z <- sel[run == s]
        ep_s <- sprintf("EP%02d", z$EP); cell_s <- sprintf("CELL%02d", z$CELL)
        ep_full <- sprintf("EP%02d", f$EP); cell_full <- sprintf("CELL%02d", f$CELL)
        r_ep <- as.numeric(ep_rows[[s]][full_factor == ep_full, ..ep_s])
        r_cell <- as.numeric(cell_rows[[s]][full_factor == cell_full, ..cell_s])
        bm <- best[mode == "EP" & subset == s & full_factor == ep_full]
        cat(sprintf("  %s: selected %s, %s | |r|(FULL %s, %s %s) = %.3f | |r|(FULL %s, %s %s) = %.3f | best match of FULL %s in %s: %s (|r| = %.3f; second %s, %.3f)\n",
                    s, ep_s, cell_s, ep_full, s, ep_s, r_ep, cell_full, s, cell_s, r_cell,
                    ep_full, s, bm$best_subset_factor, bm$best_abs_r,
                    bm$second_subset_factor, bm$second_abs_r))
    }
}

fwrite(rbindlist(ep_rows), OUT_EP, sep = "\t")
fwrite(rbindlist(cell_rows), OUT_CELL, sep = "\t")
fwrite(best, OUT_BEST, sep = "\t")
cat("\nSaved outputs with prefix GSE303006_M23_rev21b_.\n")
rev_end_log(log_con)
