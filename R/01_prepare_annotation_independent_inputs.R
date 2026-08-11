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

required_files <- c(
    atac_file,
    h3_file,
    rna_file,
    atac_peak_file,
    h3_peak_file
)

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

RUN_MACS3 <- FALSE

if (RUN_MACS3) {

    # 7列 fragment -> BED6
    system(
        paste(
            "zcat", shQuote(atac_file),
            "| awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$6}'",
            "> brain_ATAC.bed"
        )
    )

    system(
        paste(
            "zcat", shQuote(h3_file),
            "| awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$6}'",
            "> brain_H3K27ac.bed"
        )
    )

    # R を起動したshellで macs3 が見えることが前提
    system2(
        "macs3",
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
        )
    )

    system2(
        "macs3",
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
        )
    )
}

# ------------------------------------------------------------
# 3. ここが「やり直し直前」
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STOP POINT\n")
cat("M25 -> M23 のやり直しは次のファイルから開始してください:\n")
cat("02_restart_from_M23_and_after.R\n")
cat("============================================================\n")
