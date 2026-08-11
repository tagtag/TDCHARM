# ============================================================
# 02_restart_from_M23_and_after.R
# GSE303006 brain E-P analysis
#
# ★ M25 -> M23 変更後の再実行はここから
#
# 再利用するもの:
#   brain_ATAC_peaks.narrowPeak
#   brain_H3K27ac_peaks.broadPeak
#   GSE303006_brain.*.fragment.tsv.gz
#   GSE303006_brain.RNAcounts.gene.total.format.tsv.gz
#   tdg20k/geo_submit/brain_20k_tdg/
#
# 新たに M23 で作り直すもの:
#   promoter
#   enhancer candidates
#   E-P candidates
#   RNA-filtered E-P
#   enhancer x cell ATAC/H3K27ac
#   common cells
#   20-kb E-P bin pairs
#   3D summaries
#
# 出力名には M23 を明示して、M25版との混同を防ぐ。
# ============================================================

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(Matrix)

# ============================================================
# 0. Input
# ============================================================

gtf_file <- "gencode.vM23.annotation.gtf.gz"

atac_file <- "GSE303006_brain.atac.fragment.tsv.gz"
h3_file   <- "GSE303006_brain.h3k27ac.fragment.tsv.gz"
rna_file  <- "GSE303006_brain.RNAcounts.gene.total.format.tsv.gz"

atac_peak_file <- "brain_ATAC_peaks.narrowPeak"
h3_peak_file   <- "brain_H3K27ac_peaks.broadPeak"

tdg_dir <- "tdg20k/geo_submit/brain_20k_tdg"

stopifnot(
    file.exists(gtf_file),
    file.exists(atac_file),
    file.exists(h3_file),
    file.exists(rna_file),
    file.exists(atac_peak_file),
    file.exists(h3_peak_file),
    dir.exists(tdg_dir)
)

canonical_chr <- c(
    paste0("chr", 1:19),
    "chrX", "chrY"
)

# ============================================================
# 1. M23 promoter: TSS +/- 1 kb
# ============================================================

gtf <- import(gtf_file)
genes <- gtf[gtf$type == "gene"]

tss_gene <- resize(
    genes,
    width = 1,
    fix = "start"
)

promoter <- promoters(
    tss_gene,
    upstream = 1000,
    downstream = 1001
)

promoter$gene_id   <- genes$gene_id
promoter$gene_name <- genes$gene_name

promoter_df <- data.frame(
    chr       = as.character(seqnames(promoter)),
    start     = start(promoter),
    end       = end(promoter),
    gene_id   = promoter$gene_id,
    gene_name = promoter$gene_name,
    strand    = as.character(strand(promoter))
)

bed <- data.frame(
    chr    = as.character(seqnames(promoter)),
    start  = start(promoter) - 1,   # BED: 0-based start
    end    = end(promoter),
    name   = promoter$gene_name,
    score  = 0,
    strand = as.character(strand(promoter))
)

write.table(
    bed,
    file = "mm10_M23_promoter_TSS_1kb.bed",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

write.table(
    promoter_df,
    file = "mm10_M23_promoter_TSS_1kb.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("M23 genes/promoters :", length(promoter), "\n")

# ============================================================
# 2. ATAC peak intersect H3K27ac broadPeak - M23 promoter
#    MACS3 peak filesは再利用
# ============================================================

atac <- import(atac_peak_file)
h3   <- import(h3_peak_file)
prom <- import("mm10_M23_promoter_TSS_1kb.bed")

atac <- keepSeqlevels(
    atac,
    intersect(seqlevels(atac), canonical_chr),
    pruning.mode = "coarse"
)

h3 <- keepSeqlevels(
    h3,
    intersect(seqlevels(h3), canonical_chr),
    pruning.mode = "coarse"
)

prom <- keepSeqlevels(
    prom,
    intersect(seqlevels(prom), canonical_chr),
    pruning.mode = "coarse"
)

cat("ATAC peaks       :", length(atac), "\n")
cat("H3K27ac peaks    :", length(h3), "\n")
cat("M23 promoters    :", length(prom), "\n")

hit_h3 <- findOverlaps(
    atac,
    h3,
    ignore.strand = TRUE
)

active <- atac[unique(queryHits(hit_h3))]

hit_prom <- findOverlaps(
    active,
    prom,
    ignore.strand = TRUE
)

promoter_overlap <- unique(queryHits(hit_prom))

if (length(promoter_overlap) > 0) {
    enhancer <- active[-promoter_overlap]
} else {
    enhancer <- active
}

enhancer <- sort(enhancer)
enhancer$enhancer_id <- sprintf("E%06d", seq_along(enhancer))

cat("ATAC + H3K27ac              :", length(active), "\n")
cat("Removed by M23 promoter     :", length(promoter_overlap), "\n")
cat("M23 enhancer candidates     :", length(enhancer), "\n")

export(
    enhancer,
    "GSE303006_M23_brain_enhancer_candidates.bed",
    format = "BED"
)

enhancer_df <- data.frame(
    enhancer_id = enhancer$enhancer_id,
    chr         = as.character(seqnames(enhancer)),
    start       = start(enhancer),
    end         = end(enhancer),
    width       = width(enhancer)
)

fwrite(
    enhancer_df,
    "GSE303006_M23_brain_enhancer_candidates.tsv",
    sep = "\t"
)

# ============================================================
# 3. M23 enhancer x promoter: +/-1 Mb candidate E-P
# ============================================================

enhancer <- import("GSE303006_M23_brain_enhancer_candidates.bed")
enhancer$enhancer_id <- sprintf("E%06d", seq_along(enhancer))

prom <- fread("mm10_M23_promoter_TSS_1kb.tsv")
prom[, tss := floor((start + end) / 2)]

enhancer_center_gr <- resize(
    enhancer,
    width = 1,
    fix = "center"
)

enhancer_center_gr$enhancer_id <- enhancer$enhancer_id
mcols(enhancer_center_gr)$enhancer_start <- start(enhancer)
mcols(enhancer_center_gr)$enhancer_end   <- end(enhancer)

tss_gr <- GRanges(
    seqnames = prom$chr,
    ranges = IRanges(
        start = prom$tss,
        width = 1
    ),
    strand = prom$strand
)

tss_gr$gene_id   <- prom$gene_id
tss_gr$gene_name <- prom$gene_name

enhancer_center_gr <- keepSeqlevels(
    enhancer_center_gr,
    intersect(seqlevels(enhancer_center_gr), canonical_chr),
    pruning.mode = "coarse"
)

tss_gr <- keepSeqlevels(
    tss_gr,
    intersect(seqlevels(tss_gr), canonical_chr),
    pruning.mode = "coarse"
)

tss_window <- resize(
    tss_gr,
    width = 2000001,
    fix = "center"
)

hits <- findOverlaps(
    enhancer_center_gr,
    tss_window,
    ignore.strand = TRUE
)

ei <- queryHits(hits)
gi <- subjectHits(hits)

ep <- data.table(
    enhancer_id = enhancer_center_gr$enhancer_id[ei],
    enhancer_chr = as.character(seqnames(enhancer_center_gr))[ei],
    enhancer_start = enhancer_center_gr$enhancer_start[ei],
    enhancer_end = enhancer_center_gr$enhancer_end[ei],
    enhancer_center = start(enhancer_center_gr)[ei],
    gene_id = tss_gr$gene_id[gi],
    gene_name = tss_gr$gene_name[gi],
    gene_strand = as.character(strand(tss_gr))[gi],
    tss = start(tss_gr)[gi]
)

ep[, distance := enhancer_center - tss]
ep[, abs_distance := abs(distance)]

ep[, oriented_distance :=
    ifelse(
        gene_strand == "+",
        enhancer_center - tss,
        tss - enhancer_center
    )
]

ep <- ep[abs_distance <= 1000000]

ep[, pair_id := sprintf("EP%07d", seq_len(.N))]

setcolorder(
    ep,
    c(
        "pair_id",
        "enhancer_id",
        "enhancer_chr",
        "enhancer_start",
        "enhancer_end",
        "enhancer_center",
        "gene_id",
        "gene_name",
        "gene_strand",
        "tss",
        "distance",
        "abs_distance",
        "oriented_distance"
    )
)

fwrite(
    ep,
    "GSE303006_M23_brain_candidate_EP_1Mb.tsv.gz",
    sep = "\t"
)

cat("Enhancers       :", length(enhancer_center_gr), "\n")
cat("Genes/TSS       :", length(tss_gr), "\n")
cat("Candidate pairs :", nrow(ep), "\n")

# ============================================================
# 4. RNA filter: >=1% cellsで検出されたgene
# ============================================================

con <- gzfile(rna_file, "rt")
header <- strsplit(readLines(con, n = 1), "\t")[[1]]
close(con)

cells <- header[-1]
ncell <- length(cells)
min_cells <- ceiling(ncell * 0.01)

cat("RNA cells              :", ncell, "\n")
cat("1% detection threshold :", min_cells, "\n")

cmd <- sprintf(
    paste0(
        "zcat %s | ",
        "awk 'BEGIN{OFS=\"\\t\"} ",
        "NR>1 {",
        "s=0; n=0; ",
        "for(i=2;i<=NF;i++){s+=$i; if($i>0)n++} ",
        "print $1,s,n",
        "}'"
    ),
    shQuote(rna_file)
)

rna_stat <- fread(
    cmd = cmd,
    col.names = c(
        "gene_name",
        "total_count",
        "detected_cells"
    )
)

expressed_genes <- rna_stat[
    detected_cells >= min_cells
]

ep_rna <- ep[
    gene_name %in% expressed_genes$gene_name
]

# RNA matrixはgene symbolなので、
# gene_name -> gene_id が複数対応する曖昧symbolを除外
gene_name_map <- unique(
    ep_rna[, .(gene_name, gene_id)]
)

ambiguous_gene_names <- gene_name_map[
    ,
    .N,
    by = gene_name
][N > 1, gene_name]

ep_rna_unique <- ep_rna[
    !gene_name %in% ambiguous_gene_names
]

cat("All RNA genes             :", nrow(rna_stat), "\n")
cat("Expressed RNA genes       :", nrow(expressed_genes), "\n")
cat("Ambiguous gene symbols    :", length(ambiguous_gene_names), "\n")
cat("E-P before RNA filter     :", nrow(ep), "\n")
cat("E-P after RNA/symbol filter:", nrow(ep_rna_unique), "\n")
cat("Unique E-P genes          :", uniqueN(ep_rna_unique$gene_name), "\n")
cat("Unique E-P enhancers      :", uniqueN(ep_rna_unique$enhancer_id), "\n")

dup_check <- ep_rna_unique[
    ,
    .N,
    by = .(enhancer_id, gene_name)
][N > 1]

if (nrow(dup_check) != 0) {
    stop("Duplicate enhancer_id x gene_name remains after filtering.")
}

fwrite(
    ep_rna_unique,
    "GSE303006_M23_brain_candidate_EP_1Mb_RNAfiltered_unique.tsv.gz",
    sep = "\t"
)

# ============================================================
# 5. enhancer x cell sparse matrixを作る関数
# ============================================================

enh_df <- fread(
    "GSE303006_M23_brain_enhancer_candidates.tsv"
)

enh_gr <- GRanges(
    seqnames = enh_df$chr,
    ranges = IRanges(
        start = enh_df$start,
        end = enh_df$end
    )
)

enh_gr$enhancer_id <- enh_df$enhancer_id

make_enhancer_cell_matrix <- function(
    fragment_file,
    enhancer_gr,
    cells,
    chunk_size = 500000
) {

    nch <- length(enhancer_gr)
    ncell <- length(cells)

    mat <- Matrix(
        0,
        nrow = nch,
        ncol = ncell,
        sparse = TRUE
    )

    chr_use <- unique(as.character(seqnames(enhancer_gr)))
    con <- gzfile(fragment_file, "rt")
    on.exit(close(con), add = TRUE)

    chunk_no <- 0

    repeat {

        lines <- readLines(
            con,
            n = chunk_size
        )

        if (length(lines) == 0) break

        chunk_no <- chunk_no + 1

        dt <- fread(
            text = paste(lines, collapse = "\n"),
            header = FALSE,
            select = c(1, 2, 3, 7)
        )

        setnames(
            dt,
            c("chr", "start0", "end", "cell")
        )

        dt <- dt[chr %in% chr_use]
        dt[, cell_index := match(cell, cells)]
        dt <- dt[!is.na(cell_index)]

        if (nrow(dt) == 0) next

        frag_gr <- GRanges(
            seqnames = dt$chr,
            ranges = IRanges(
                start = dt$start0 + 1,
                end = dt$end
            )
        )

        hits <- findOverlaps(
            enhancer_gr,
            frag_gr,
            ignore.strand = TRUE
        )

        if (length(hits) > 0) {

            ii <- queryHits(hits)
            jj <- dt$cell_index[subjectHits(hits)]

            add <- sparseMatrix(
                i = ii,
                j = jj,
                x = 1,
                dims = c(nch, ncell)
            )

            mat <- mat + add
        }

        cat(
            "chunk", chunk_no,
            "lines", length(lines),
            "nnzero", nnzero(mat),
            "\n"
        )
    }

    rownames(mat) <- enhancer_gr$enhancer_id
    colnames(mat) <- cells

    mat
}

# ============================================================
# 6. M23 enhancer x cell ATAC / H3K27ac
# ============================================================

ATAC <- make_enhancer_cell_matrix(
    atac_file,
    enhancer_gr = enh_gr,
    cells = cells
)

saveRDS(
    ATAC,
    "GSE303006_M23_brain_enhancer_by_cell_ATAC.rds"
)

H3K27ac <- make_enhancer_cell_matrix(
    h3_file,
    enhancer_gr = enh_gr,
    cells = cells
)

saveRDS(
    H3K27ac,
    "GSE303006_M23_brain_enhancer_by_cell_H3K27ac.rds"
)

cat("ATAC dim      :", paste(dim(ATAC), collapse = " x "), "\n")
cat("ATAC nnzero   :", nnzero(ATAC), "\n")
cat("H3K27ac dim   :", paste(dim(H3K27ac), collapse = " x "), "\n")
cat("H3K27ac nnzero:", nnzero(H3K27ac), "\n")

# ============================================================
# 7. ATAC + H3K27ac 共通cells
# ============================================================

atac_sum <- Matrix::colSums(ATAC)
h3_sum   <- Matrix::colSums(H3K27ac)

common_cells_AH <- cells[
    atac_sum > 0 &
    h3_sum > 0
]

qc <- data.table(
    cell = cells,
    ATAC = as.numeric(atac_sum),
    H3K27ac = as.numeric(h3_sum)
)

qc[, ATAC_zero := ATAC == 0]
qc[, H3K27ac_zero := H3K27ac == 0]
qc[, keep_AH := ATAC > 0 & H3K27ac > 0]

fwrite(
    qc,
    "GSE303006_M23_brain_cell_QC_ATAC_H3K27ac.tsv",
    sep = "\t"
)

ATAC_AH <- ATAC[
    ,
    common_cells_AH,
    drop = FALSE
]

H3K27ac_AH <- H3K27ac[
    ,
    common_cells_AH,
    drop = FALSE
]

saveRDS(
    common_cells_AH,
    "GSE303006_M23_brain_common_cells_ATAC_H3K27ac.rds"
)

saveRDS(
    ATAC_AH,
    "GSE303006_M23_brain_ATAC_enhancer_by_cell_AHfiltered.rds"
)

saveRDS(
    H3K27ac_AH,
    "GSE303006_M23_brain_H3K27ac_enhancer_by_cell_AHfiltered.rds"
)

save(
    ATAC_AH,
    H3K27ac_AH,
    common_cells_AH,
    ep_rna_unique,
    file = "GSE303006_M23_brain_EP_analysis_stage_AH.RData"
)

cat("ATAC/H3K27ac common cells :", length(common_cells_AH), "\n")

# ============================================================
# 8. E-P -> 20 kb bin pair
# ============================================================

ep3 <- copy(ep_rna_unique)

ep3[, enhancer_bin :=
    floor(enhancer_center / 20000) * 20000
]

ep3[, promoter_bin :=
    floor(tss / 20000) * 20000
]

ep3[, same_20k_bin :=
    enhancer_bin == promoter_bin
]

ep_bins <- unique(
    ep3[, .(
        enhancer_chr,
        enhancer_bin,
        promoter_bin,
        same_20k_bin
    )]
)

ep_bins[, binpair_id := .I]

ep3[
    ep_bins,
    on = .(
        enhancer_chr,
        enhancer_bin,
        promoter_bin
    ),
    binpair_id := i.binpair_id
]

cat("E-P pairs            :", nrow(ep3), "\n")
cat("Unique 20kb binpairs :", nrow(ep_bins), "\n")
cat("Same-bin E-P         :", sum(ep3$same_20k_bin), "\n")
cat("Same-bin binpairs    :", sum(ep_bins$same_20k_bin), "\n")

saveRDS(
    ep3,
    "GSE303006_M23_brain_EP_with_20kb_binpair.rds"
)

saveRDS(
    ep_bins,
    "GSE303006_M23_brain_unique_20kb_binpairs.rds"
)

# ============================================================
# 9. 3DファイルとATAC/H3K27ac共通cell
# ============================================================

tdg_files <- list.files(
    tdg_dir,
    pattern = "\\.20k\\.3dg\\.gz$",
    full.names = TRUE
)

tdg_cells <- sub(
    "\\.20k\\.3dg\\.gz$",
    "",
    basename(tdg_files)
)

file_by_cell <- setNames(
    tdg_files,
    tdg_cells
)

common_cells_AH3D <- common_cells_AH[
    common_cells_AH %in% tdg_cells
]

cat("ATAC/H3K27ac cells :", length(common_cells_AH), "\n")
cat("Also having 3D     :", length(common_cells_AH3D), "\n")

missing_3D_cells <- setdiff(
    common_cells_AH,
    tdg_cells
)

if (length(missing_3D_cells) > 0) {
    cat(
        "Cells without 3D:\n",
        paste(missing_3D_cells, collapse = " "),
        "\n"
    )
}

saveRDS(
    common_cells_AH3D,
    "GSE303006_M23_brain_common_cells_ATAC_H3K27ac_3D.rds"
)

ATAC_AH3D <- ATAC_AH[
    ,
    common_cells_AH3D,
    drop = FALSE
]

H3K27ac_AH3D <- H3K27ac_AH[
    ,
    common_cells_AH3D,
    drop = FALSE
]

saveRDS(
    ATAC_AH3D,
    "GSE303006_M23_brain_ATAC_enhancer_by_cell_AH3D.rds"
)

saveRDS(
    H3K27ac_AH3D,
    "GSE303006_M23_brain_H3K27ac_enhancer_by_cell_AH3D.rds"
)

# ============================================================
# 10. 3Dで使える異なる20 kb binだけに限定
# ============================================================

ep_bins_3D <- ep_bins[
    same_20k_bin == FALSE
]

ep3_3D <- ep3[
    same_20k_bin == FALSE
]

cat("3D bin-pairs :", nrow(ep_bins_3D), "\n")
cat("3D E-P pairs :", nrow(ep3_3D), "\n")

# ============================================================
# 11. 必要な20 kb lociを一意化
# ============================================================

loci <- unique(
    rbind(
        ep_bins_3D[, .(
            chr = enhancer_chr,
            pos = enhancer_bin
        )],
        ep_bins_3D[, .(
            chr = enhancer_chr,
            pos = promoter_bin
        )]
    )
)

loci[, locus_id := .I]
loci[, locus_key := paste(chr, pos, sep = ":")]
loci[, key_mat := paste(chr, pos, "mat", sep = ":")]
loci[, key_pat := paste(chr, pos, "pat", sep = ":")]

locus_index <- setNames(
    loci$locus_id,
    loci$locus_key
)

ep_bins_3D[, e_locus :=
    unname(
        locus_index[
            paste(enhancer_chr, enhancer_bin, sep = ":")
        ]
    )
]

ep_bins_3D[, p_locus :=
    unname(
        locus_index[
            paste(enhancer_chr, promoter_bin, sep = ":")
        ]
    )
]

stopifnot(
    !anyNA(ep_bins_3D$e_locus),
    !anyNA(ep_bins_3D$p_locus)
)

cat("Unique 20kb loci :", nrow(loci), "\n")

saveRDS(
    loci,
    "GSE303006_M23_brain_20kb_loci_for_EP.rds"
)

saveRDS(
    ep_bins_3D,
    "GSE303006_M23_brain_20kb_binpairs_indexed.rds"
)

# ============================================================
# 12. 3dg reader
# ============================================================

read_3dg <- function(file) {

    d <- fread(
        cmd = paste("zcat", shQuote(file)),
        header = FALSE,
        col.names = c(
            "chrom_phase",
            "position",
            "x",
            "y",
            "z"
        )
    )

    d[, chr :=
        sub("\\((mat|pat)\\)$", "", chrom_phase)
    ]

    d[, phase :=
        sub("^.*\\((mat|pat)\\)$", "\\1", chrom_phase)
    ]

    d[, key :=
        paste(chr, position, phase, sep = ":")
    ]

    d
}

# ============================================================
# 13. 高速3D距離関数
#     cis E-Pなので mat-mat / pat-pat のみを使い、
#     unphased解析用に両haplotypeの平均距離を使う。
# ============================================================

calc_3D_fast <- function(bp, loci, tdg) {

    idx_mat <- match(
        loci$key_mat,
        tdg$key
    )

    idx_pat <- match(
        loci$key_pat,
        tdg$key
    )

    em <- idx_mat[bp$e_locus]
    pm <- idx_mat[bp$p_locus]

    epat <- idx_pat[bp$e_locus]
    ppat <- idx_pat[bp$p_locus]

    n <- nrow(bp)

    d_mat <- rep(NA_real_, n)
    d_pat <- rep(NA_real_, n)

    ok_mat <- !is.na(em) & !is.na(pm)

    d_mat[ok_mat] <- sqrt(
        (tdg$x[em[ok_mat]] - tdg$x[pm[ok_mat]])^2 +
        (tdg$y[em[ok_mat]] - tdg$y[pm[ok_mat]])^2 +
        (tdg$z[em[ok_mat]] - tdg$z[pm[ok_mat]])^2
    )

    ok_pat <- !is.na(epat) & !is.na(ppat)

    d_pat[ok_pat] <- sqrt(
        (tdg$x[epat[ok_pat]] - tdg$x[ppat[ok_pat]])^2 +
        (tdg$y[epat[ok_pat]] - tdg$y[ppat[ok_pat]])^2 +
        (tdg$z[epat[ok_pat]] - tdg$z[ppat[ok_pat]])^2
    )

    d_mean <- rowMeans(
        cbind(d_mat, d_pat),
        na.rm = TRUE
    )

    d_mean[is.nan(d_mean)] <- NA_real_

    d_mean
}

# ============================================================
# 14. 1 cellで3D計算をテスト
# ============================================================

test_cell <- common_cells_AH3D[1]

tdg_test <- read_3dg(
    file_by_cell[[test_cell]]
)

v_test <- calc_3D_fast(
    ep_bins_3D,
    loci,
    tdg_test
)

cat("Test 3D cell     :", test_cell, "\n")
cat("Test 3D coverage :", mean(!is.na(v_test)), "\n")
print(summary(v_test))

rm(tdg_test, v_test)
gc()

# ============================================================
# 15. 全cellについて3D summary
#     checkpoint付き。M23専用checkpoint名。
# ============================================================

cells_run <- common_cells_AH3D

stopifnot(
    all(cells_run %in% names(file_by_cell))
)

checkpoint_file <- "GSE303006_M23_3D_summary_checkpoint.rds"
n_bp <- nrow(ep_bins_3D)

if (file.exists(checkpoint_file)) {

    chk <- readRDS(checkpoint_file)

    # M23の現在のbinpair数と一致しないcheckpointは使わない
    if (length(chk$n_obs) != n_bp) {
        stop(
            "Checkpoint dimensions do not match current M23 binpairs. ",
            "Delete/rename ", checkpoint_file, " and restart."
        )
    }

    n_obs   <- chk$n_obs
    sum_d   <- chk$sum_d
    sumsq_d <- chk$sumsq_d

    start_k <- chk$k + 1L

    cat(
        "Restarting from cell",
        start_k,
        "of",
        length(cells_run),
        "\n"
    )

} else {

    n_obs   <- integer(n_bp)
    sum_d   <- numeric(n_bp)
    sumsq_d <- numeric(n_bp)

    start_k <- 1L

    cat("Starting 3D summary from cell 1\n")
}

if (start_k <= length(cells_run)) {

    for (k in seq.int(
        from = start_k,
        to = length(cells_run)
    )) {

        cell <- cells_run[k]

        tdg <- read_3dg(
            file_by_cell[[cell]]
        )

        v <- calc_3D_fast(
            ep_bins_3D,
            loci,
            tdg
        )

        ok <- !is.na(v)

        n_obs[ok] <- n_obs[ok] + 1L
        sum_d[ok] <- sum_d[ok] + v[ok]
        sumsq_d[ok] <- sumsq_d[ok] + v[ok]^2

        if (k %% 25 == 0) {
            cat(
                "processed",
                k,
                "/",
                length(cells_run),
                " cell =",
                cell,
                " coverage =",
                round(mean(ok), 4),
                "\n"
            )
        }

        if (k %% 50 == 0) {

            saveRDS(
                list(
                    k = k,
                    n_obs = n_obs,
                    sum_d = sum_d,
                    sumsq_d = sumsq_d
                ),
                checkpoint_file
            )

            gc()
        }
    }
}

saveRDS(
    list(
        k = length(cells_run),
        n_obs = n_obs,
        sum_d = sum_d,
        sumsq_d = sumsq_d
    ),
    checkpoint_file
)

# ============================================================
# 16. 3D summary statistics
# ============================================================

ep_bins_3D[, n_3D := n_obs]

ep_bins_3D[, coverage_3D :=
    n_3D / length(cells_run)
]

ep_bins_3D[, mean_distance_3D :=
    fifelse(
        n_3D > 0,
        sum_d / n_3D,
        NA_real_
    )
]

ep_bins_3D[, sd_distance_3D :=
    fifelse(
        n_3D > 1,
        sqrt(
            pmax(
                0,
                (sumsq_d - sum_d^2 / n_3D) /
                (n_3D - 1)
            )
        ),
        NA_real_
    )
]

cat("\n3D coverage summary\n")
print(summary(ep_bins_3D$coverage_3D))

cat("\nMean 3D distance summary\n")
print(summary(ep_bins_3D$mean_distance_3D))

cat("\nSD 3D distance summary\n")
print(summary(ep_bins_3D$sd_distance_3D))

cat(
    "\nZero-coverage binpairs :",
    sum(ep_bins_3D$coverage_3D == 0),
    "\n"
)

saveRDS(
    ep_bins_3D,
    "GSE303006_M23_brain_20kb_3D_summary_complete.rds"
)

fwrite(
    ep_bins_3D,
    "GSE303006_M23_brain_20kb_3D_summary_complete.tsv.gz",
    sep = "\t"
)

cat("\n============================================================\n")
cat("M23 rebuild through 3D summary completed.\n")
cat("============================================================\n")
