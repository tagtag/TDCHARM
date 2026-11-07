# ============================================================
# tests/make_synthetic_workdir.R
#
# Builds a SMALL synthetic working directory with the same file
# names and object structures as the real GSE303006 intermediate
# files expected by scripts 05-11 and the revision scripts 19-27.
# Used only to check that the scripts run end to end and that the
# revision code reproduces the algebra of the original pipeline.
# The numbers produced from it have no biological meaning.
#
# Usage:
#   Rscript tests/make_synthetic_workdir.R /path/to/empty_dir
# ============================================================

suppressPackageStartupMessages({
    library(data.table); library(Matrix); library(rhdf5)
})

args <- commandArgs(trailingOnly = TRUE)
wd <- if (length(args)) args[1] else "synthetic_work"
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
set.seed(20261103)

types <- c("Inh Ndnf/Lamp5", "Ex L3/4 IT", "Ex L2/3 IT", "Astro")
n_per_rep <- 70L
cells <- unlist(lapply(c("R1", "R2", "R3"), function(r)
    sprintf("%s_cell%03d", r, seq_len(n_per_rep))))
nC <- length(cells)
ctype <- sample(types, nC, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))

n_enh_all <- 320L; n_enh <- 260L; n_gene <- 140L
enh_all <- sprintf("E%06d", seq_len(n_enh_all))
enh <- enh_all[seq_len(n_enh)]
genes <- sprintf("Gene%03d", seq_len(n_gene))

# E-P pairs: each final enhancer linked to 1-4 genes
ep <- rbindlist(lapply(seq_len(n_enh), function(i) {
    g <- sample(genes, sample(1:4, 1))
    data.table(enhancer_id = enh[i], gene_name = g)
}))
ep[, enhancer_chr := "chr1"]
ep[, enhancer_start := match(enhancer_id, enh) * 5000L]
ep[, enhancer_end := enhancer_start + 500L]
ep[, enhancer_center := enhancer_start + 250L]
ep[, tss := match(gene_name, genes) * 9000L]
ep[, distance := enhancer_center - tss]
ep[, abs_distance := abs(distance)]
ep[, pair_id := sprintf("EP%07d", seq_len(.N))]
ep[, binpair_id := as.integer(factor(paste(enhancer_start %/% 20000L,
                                           tss %/% 20000L)))]
saveRDS(ep, "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds")
bp <- sort(unique(ep$binpair_id))

# Latent programs: program 1 separates Ndnf/Lamp5 from L3/4 IT
prog <- cbind(
    p1 = ifelse(ctype == types[1], 1.5, ifelse(ctype == types[2], -1.5, 0)),
    p2 = ifelse(ctype == types[4], 1.2, -0.3),
    p3 = rnorm(nC))
rep_shift <- c(R1 = 0, R2 = 0.15, R3 = -0.1)[substr(cells, 1, 2)]

counts <- function(nf, load_scale = 1, mu = 1,
                   L = matrix(rnorm(nf * 3, sd = load_scale), nf, 3)) {
    eta <- log(mu) + L %*% t(prog) * 0.6 + rep_shift
    m <- matrix(rpois(nf * nC, exp(eta)), nf, nC)
    Matrix(m, sparse = TRUE)
}
# ATAC and H3K27ac share enhancer loadings (correlated modalities)
L_enh <- matrix(rnorm(n_enh_all * 3), n_enh_all, 3)
ATAC <- counts(n_enh_all, mu = 0.8, L = L_enh); dimnames(ATAC) <- list(enh_all, cells)
H3 <- counts(n_enh_all, mu = 0.6, L = L_enh); dimnames(H3) <- list(enh_all, cells)
saveRDS(ATAC, "GSE303006_M23_brain_ATAC_enhancer_by_cell_AH3D.rds")
saveRDS(H3, "GSE303006_M23_brain_H3K27ac_enhancer_by_cell_AH3D.rds")

RNAc <- counts(n_gene, mu = 2); dimnames(RNAc) <- list(genes, cells)
lib <- colSums(RNAc)
RNA <- RNAc %*% Diagonal(x = 1e4 / lib); RNA <- as(RNA, "dgCMatrix")
RNA@x <- log1p(RNA@x); dimnames(RNA) <- list(genes, cells)
saveRDS(RNA, "GSE303006_M23_brain_RNA_logCP10K_gene_by_cell.rds")

qc <- data.table(cell = cells, ATAC_count = colSums(ATAC[enh, ]),
                 H3K27ac_count = colSums(H3[enh, ]),
                 RNA_library_size = lib,
                 RNA_detected_genes = colSums(RNAc > 0), keep_final = TRUE)
fwrite(qc, "GSE303006_M23_brain_final_cell_QC.tsv", sep = "\t")
saveRDS(cells, "GSE303006_M23_brain_final_4258_cells.rds")

# 3D distances (bin pair x cell), some missing
nb <- length(bp)
Lb <- matrix(rnorm(nb * 3, sd = 0.3), nb, 3)
logd <- 1 + Lb %*% t(prog) + matrix(rnorm(nb * nC, sd = 0.4), nb, nC)
D <- exp(logd)
D[matrix(runif(nb * nC) < 0.05, nb, nC)] <- NA
h5 <- "GSE303006_M23_brain_3D_binpair_by_cell.h5"
if (file.exists(h5)) file.remove(h5)
h5createFile(h5)
h5write(as.integer(bp), h5, "binpair_id")
h5write(cells, h5, "cell")
h5createDataset(h5, "distance", dims = c(nb, nC), storage.mode = "double",
                chunk = c(min(nb, 50L), nC))
h5write(D, h5, "distance")
h5closeAll()

meta <- data.table(cellname = cells, celltype = ctype,
                   sample = substr(cells, 1, 2))
fwrite(meta, "GSE303006_charm_metadata_qcpass.tsv.gz", sep = "\t")

cat("Synthetic working directory written to", normalizePath("."), "\n")
cat("cells:", nC, " enhancers:", n_enh, " genes:", n_gene,
    " E-P pairs:", nrow(ep), " bin pairs:", nb, "\n")
