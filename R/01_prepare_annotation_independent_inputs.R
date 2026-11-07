# ============================================================
# 01_before_M23_restart.R
# GSE303006 brain E-P analysis
#
# このファイルは GENCODE M23/M25 に依存しない部分です。
# 現在すでに以下が存在するなら、通常は再実行不要です。
#
#   brain_ATAC_peaks.narrowPeak
#   brain_H3K27ac_peaks.broadPeak
#   GSE303006_brain.RNAcounts.gene.total.format.tsv.gz
#   tdg20k/geo_submit/brain_20k_tdg/
#
# M25 -> M23 の変更でやり直すのは 02_restart_from_M23_and_after.R
# からです。
# ============================================================

library(data.table)
library(GenomicRanges)
library(rtracklayer)

atac_file <- "GSE303006_brain.atac.fragment.tsv.gz"
h3_file   <- "GSE303006_brain.h3k27ac.fragment.tsv.gz"
rna_file  <- "GSE303006_brain.RNAcounts.gene.total.format.tsv.gz"

atac_peak_file <- "brain_ATAC_peaks.narrowPeak"
h3_peak_file   <- "brain_H3K27ac_peaks.broadPeak"

tdg_dir <- "tdg20k/geo_submit/brain_20k_tdg"

# ------------------------------------------------------------
# 1. 現在の再利用可能ファイルを確認
# ------------------------------------------------------------

# Revision: MACS3 peak calling is run when the peak files are absent
# (full recomputation from the GEO files) or when TDCHARM_RUN_MACS3=1.
# Otherwise the existing peak files are reused, as before.
RUN_MACS3 <- identical(Sys.getenv("TDCHARM_RUN_MACS3"), "1") ||
    !file.exists(atac_peak_file) || !file.exists(h3_peak_file)

required_files <- c(
    atac_file,
    h3_file,
    rna_file
)
if (!RUN_MACS3) {
    required_files <- c(required_files, atac_peak_file, h3_peak_file)
}

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    stop(
        "Missing reusable files:\n",
        paste(missing_files, collapse = "\n")
    )
}

if (!dir.exists(tdg_dir)) {
    stop("3D directory not found: ", tdg_dir)
}

cat("Annotation-independent inputs are present.\n")
cat("ATAC peak file    :", atac_peak_file, "\n")
cat("H3K27ac peak file :", h3_peak_file, "\n")
cat("3D directory      :", tdg_dir, "\n")

# ------------------------------------------------------------
# 2. 参考：MACS3 peak を作った手順
#    M23/M25 annotation には依存しないので、既存peakを再利用する。
#    再作成する必要がある場合だけ RUN_MACS3 <- TRUE にする。
# ------------------------------------------------------------

if (RUN_MACS3) {

    # MACS3 executable: TDCHARM_MACS3 (full path, e.g. from a separate
    # conda environment) or "macs3" on PATH.
    MACS3_BIN <- Sys.getenv("TDCHARM_MACS3", unset = "macs3")
    if (!nzchar(Sys.which(MACS3_BIN)) && !file.exists(MACS3_BIN)) {
        stop("Peak files are missing and MACS3 was not found (", MACS3_BIN, "). ",
             "Install it, e.g. in its own conda environment:\n",
             "  conda create -n macs3 -c conda-forge -c bioconda python=3.12 macs3=3.1.0\n",
             "and set TDCHARM_MACS3=<conda base>/envs/macs3/bin/macs3.")
    }
    # Every command, the MACS3 version and the resulting peak counts are
    # written to MACS3_LOG; the MACS3 console output of each call goes to
    # <name>_macs3_callpeak.log. Both are kept with the outputs.
    MACS3_LOG <- "GSE303006_M23_MACS3_peak_calling_log.txt"
    macs3_version <- paste(system2(MACS3_BIN, "--version", stdout = TRUE, stderr = TRUE),
                           collapse = " ")
    writeLines(c(
        "MACS3 peak calling (01_prepare_annotation_independent_inputs.R)",
        paste0("Run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("Working directory: ", getwd()),
        paste0("MACS3 executable: ", normalizePath(Sys.which(MACS3_BIN), mustWork = FALSE)),
        paste0("MACS3 version: ", macs3_version),
        ""), MACS3_LOG)
    log_line <- function(...) cat(..., "\n", sep = "", file = MACS3_LOG, append = TRUE)
    cat("Running MACS3 peak calling (", MACS3_BIN, ": ", macs3_version, ")\n", sep = "")

    run_or_stop <- function(cmd, args = character(), log_file = NULL) {
        shown <- if (length(args)) paste(c(cmd, args), collapse = " ") else cmd
        log_line("$ ", shown)
        t0 <- Sys.time()
        st <- if (length(args)) {
            if (is.null(log_file)) system2(cmd, args)
            else system2(cmd, args, stdout = log_file, stderr = log_file)
        } else {
            system(cmd)
        }
        log_line("  exit status ", st, ", ",
                 format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
        if (!identical(as.integer(st), 0L)) stop("Command failed (status ", st, "): ", shown)
        invisible(NULL)
    }

    # 7列 fragment -> BED6
    run_or_stop(
        paste(
            "zcat", shQuote(atac_file),
            "| awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$6}'",
            "> brain_ATAC.bed"
        )
    )

    run_or_stop(
        paste(
            "zcat", shQuote(h3_file),
            "| awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$6}'",
            "> brain_H3K27ac.bed"
        )
    )

    run_or_stop(
        MACS3_BIN,
        args = c(
            "callpeak",
            "-t", "brain_ATAC.bed",
            "-f", "BED",
            "-g", "mm",
            "-n", "brain_ATAC",
            "--nomodel",
            "--shift", "-50",
            "--extsize", "100",
            "--keep-dup", "all",
            "-q", "0.01"
        ),
        log_file = "brain_ATAC_macs3_callpeak.log"
    )

    run_or_stop(
        MACS3_BIN,
        args = c(
            "callpeak",
            "-t", "brain_H3K27ac.bed",
            "-f", "BED",
            "-g", "mm",
            "-n", "brain_H3K27ac",
            "--nomodel",
            "--shift", "-50",
            "--extsize", "100",
            "--keep-dup", "all",
            "--broad",
            "--broad-cutoff", "0.05"
        ),
        log_file = "brain_H3K27ac_macs3_callpeak.log"
    )

    stopifnot(file.exists(atac_peak_file), file.exists(h3_peak_file))
    n_atac <- length(readLines(atac_peak_file))
    n_h3 <- length(readLines(h3_peak_file))
    log_line("")
    log_line("Peaks: ", atac_peak_file, " = ", n_atac, ", ", h3_peak_file, " = ", n_h3)
    cat("MACS3 peak files written:", atac_peak_file, "(", n_atac, "peaks ),",
        h3_peak_file, "(", n_h3, "peaks )\n")
    cat("Log: ", MACS3_LOG, "\n", sep = "")
}

# ------------------------------------------------------------
# 3. ここが「やり直し直前」
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STOP POINT\n")
cat("M25 -> M23 のやり直しは次のファイルから開始してください:\n")
cat("02_restart_from_M23_and_after.R\n")
cat("============================================================\n")
