
# ============================================================
# 05_stage1_independentK20.R
# GSE303006 / GENCODE M23
#
# Stage 1:
#   Each modality -> K=20 components independently.
#
# Stage-2 conceptual tensor:
#   E-P x cell x K x modality
#   730,969 x 4,258 x 20 x 4
#
# IMPORTANT:
#   k=1..20 means the singular-value rank WITHIN each modality.
#   Modality remains a separate 4th mode.
#
# Preprocessing:
#   ATAC/H3K27ac: logCP10K -> feature-wise z
#   RNA: existing logCP10K -> feature-wise z
#   3D: -log(distance) -> binpair-wise z -> missing z = 0
#
# Equal modality weighting:
#   The retained K-component reconstruction of each modality
#   is scaled to Frobenius norm 1.
# ============================================================

library(data.table)
library(Matrix)
library(rhdf5)
library(RSpectra)

K <- 20L
SCALE_FACTOR <- 10000
EPS_SD <- 1e-8

# randomized SVD settings for 3D
OVERSAMPLE <- 10L
N_POWER <- 1L
BLOCK_3D <- 2000L
RANDOM_SEED <- 20260811L

EP_FILE <- "GSE303006_M23_brain_EP_final_3Dcov90_enhQC005.rds"
CELL_FILE <- "GSE303006_M23_brain_final_4258_cells.rds"
QC_FILE <- "GSE303006_M23_brain_final_cell_QC.tsv"

RNA_FILE <- "GSE303006_M23_brain_RNA_logCP10K_gene_by_cell.rds"
ATAC_FILE <- "GSE303006_M23_brain_ATAC_enhancer_by_cell_AH3D.rds"
H3_FILE <- "GSE303006_M23_brain_H3K27ac_enhancer_by_cell_AH3D.rds"
H5_3D <- "GSE303006_M23_brain_3D_binpair_by_cell.h5"

# ---------- final cells ----------
if (!file.exists(CELL_FILE)) {
    if (!file.exists(QC_FILE))
        stop("Neither final cell RDS nor QC TSV exists.")
    qc <- fread(QC_FILE)
    if (!all(c("cell","keep_final") %in% names(qc)))
        stop("QC TSV lacks cell/keep_final.")
    final_cells <- qc[keep_final == TRUE, as.character(cell)]
    saveRDS(final_cells, CELL_FILE)
} else {
    final_cells <- as.character(readRDS(CELL_FILE))
}

# ---------- E-P mapping ----------
ep <- as.data.table(readRDS(EP_FILE))

enh <- sort(unique(as.character(ep$enhancer_id)))
genes <- sort(unique(as.character(ep$gene_name)))
bp <- sort(unique(as.integer(ep$binpair_id)))

cat("E-P pairs :", nrow(ep), "\n")
cat("Enhancers :", length(enh), "\n")
cat("Genes     :", length(genes), "\n")
cat("Bin-pairs :", length(bp), "\n")
cat("Cells     :", length(final_cells), "\n")

# ---------- load 3 in-memory modalities ----------
ATAC_all <- readRDS(ATAC_FILE)
H3_all <- readRDS(H3_FILE)
RNA_all <- readRDS(RNA_FILE)

stopifnot(all(final_cells %in% colnames(ATAC_all)))
stopifnot(all(final_cells %in% colnames(H3_all)))
stopifnot(all(final_cells %in% colnames(RNA_all)))
stopifnot(all(enh %in% rownames(ATAC_all)))
stopifnot(all(enh %in% rownames(H3_all)))
stopifnot(all(genes %in% rownames(RNA_all)))

# Depth uses all enhancer-candidate rows, not only final 42,669.
ATAC_depth <- as.numeric(Matrix::colSums(ATAC_all[, final_cells, drop=FALSE]))
H3_depth <- as.numeric(Matrix::colSums(H3_all[, final_cells, drop=FALSE]))
names(ATAC_depth) <- final_cells
names(H3_depth) <- final_cells

ATAC <- ATAC_all[enh, final_cells, drop=FALSE]
H3 <- H3_all[enh, final_cells, drop=FALSE]
RNA <- RNA_all[genes, final_cells, drop=FALSE]

rm(ATAC_all,H3_all,RNA_all); gc()

# ---------- logCP10K ----------
logcp10k <- function(X, depth) {
    if (any(depth <= 0)) stop("Zero depth.")
    Y <- X %*% Diagonal(x=SCALE_FACTOR/as.numeric(depth))
    Y <- as(Y, "dgCMatrix")
    Y@x <- log1p(Y@x)
    rownames(Y) <- rownames(X)
    colnames(Y) <- colnames(X)
    Y
}

ATAC <- logcp10k(ATAC, ATAC_depth)
H3 <- logcp10k(H3, H3_depth)

# ---------- feature z stats ----------
zstats <- function(X) {
    n <- ncol(X)
    mu <- as.numeric(Matrix::rowMeans(X))
    sq <- as.numeric(Matrix::rowSums(X*X))
    vr <- pmax((sq - n*mu^2)/(n-1), 0)
    sd <- sqrt(vr)
    good <- is.finite(sd) & sd > EPS_SD
    list(mu=mu, sd=sd, good=good)
}

# ---------- sparse truncated SVD ----------
# X is feature x cell.
# RSpectra runs on t(X): cell x feature.
# center/scale are feature-wise and handled implicitly.
svd_sparse_z <- function(X, K, label) {

    st <- zstats(X)
    good <- st$good

    cat(label, "nonconstant features:", sum(good), "/", length(good), "\n")

    A <- Matrix::t(X[good, , drop=FALSE])

    s <- RSpectra::svds(
        A,
        k=K,
        nu=K,
        nv=K,
        opts=list(
            center=st$mu[good],
            scale=st$sd[good],
            tol=1e-8,
            maxitr=2000
        )
    )

    # t(Z) = U d V'
    # For original Z(feature x cell):
    # feature loading = V; cell loading = U.
    Vcell <- s$u
    Ufeat_good <- s$v
    d <- as.numeric(s$d)

    # deterministic sign
    for (k in seq_len(K)) {
        j <- which.max(abs(Vcell[,k]))
        if (Vcell[j,k] < 0) {
            Vcell[,k] <- -Vcell[,k]
            Ufeat_good[,k] <- -Ufeat_good[,k]
        }
    }

    # Equal weight across modalities:
    # retained K-component Frobenius norm -> 1.
    normK <- sqrt(sum(d^2))
    d_scaled <- d/normK

    B <- matrix(0, nrow(X), K)
    B[good,] <- sweep(Ufeat_good, 2, d_scaled, "*")

    rownames(B) <- rownames(X)
    colnames(B) <- sprintf("K%02d", seq_len(K))
    rownames(Vcell) <- colnames(X)
    colnames(Vcell) <- colnames(B)

    total_fro2 <- sum(good) * (ncol(X)-1)
    explained <- d^2/total_fro2

    list(
        B=B,
        V=Vcell,
        d=d,
        d_scaled=d_scaled,
        explained=explained,
        stats=st,
        normK=normK
    )
}

cat("\nATAC SVD...\n")
SVD_A <- svd_sparse_z(ATAC, K, "ATAC")

cat("\nH3K27ac SVD...\n")
SVD_H <- svd_sparse_z(H3, K, "H3K27ac")

cat("\nRNA SVD...\n")
SVD_R <- svd_sparse_z(RNA, K, "RNA")

# ============================================================
# 3D preparation
# ============================================================

h5_bp <- as.integer(h5read(H5_3D, "binpair_id"))
h5_cells <- as.character(h5read(H5_3D, "cell"))

rows3D <- match(bp, h5_bp)
cols3D <- match(final_cells, h5_cells)

stopifnot(!anyNA(rows3D), !anyNA(cols3D))

nBP <- length(bp)
nC <- length(final_cells)

blocks <- split(
    seq_len(nBP),
    ceiling(seq_len(nBP)/BLOCK_3D)
)

Dmu <- numeric(nBP)
Dsd <- numeric(nBP)
Dn <- integer(nBP)
Dgood <- logical(nBP)

# Read raw 3D -> standardized block.
get_z3d <- function(ii) {

    d <- as.matrix(
        h5read(
            H5_3D,
            "distance",
            index=list(rows3D[ii], cols3D)
        )
    )

    d[!is.finite(d) | d <= 0] <- NA_real_
    q <- -log(d)

    cen <- sweep(q, 1, Dmu[ii], "-")
    z <- matrix(0, length(ii), nC)

    good <- Dgood[ii]

    if (any(good)) {
        z[good,] <- sweep(
            cen[good,,drop=FALSE],
            1,
            Dsd[ii[good]],
            "/"
        )
        z[!is.finite(z)] <- 0
    }

    z
}

# ---------- first 3D pass: mean / sd ----------
cat("\n3D z-score statistics...\n")

for (ib in seq_along(blocks)) {

    ii <- blocks[[ib]]

    d <- as.matrix(
        h5read(
            H5_3D,
            "distance",
            index=list(rows3D[ii], cols3D)
        )
    )

    d[!is.finite(d) | d <= 0] <- NA_real_
    q <- -log(d)

    nobs <- rowSums(!is.na(q))
    mu <- rowMeans(q, na.rm=TRUE)
    mu[!is.finite(mu)] <- 0

    cen <- sweep(q,1,mu,"-")
    ss <- rowSums(cen^2, na.rm=TRUE)

    sd <- numeric(length(ii))
    ok <- nobs > 1
    sd[ok] <- sqrt(ss[ok]/(nobs[ok]-1))

    good <- is.finite(sd) & sd > EPS_SD

    Dmu[ii] <- mu
    Dsd[ii] <- sd
    Dn[ii] <- nobs
    Dgood[ii] <- good

    if (ib%%10==0 || ib==length(blocks))
        cat("3D stats",ib,"/",length(blocks),"\n")

    rm(d,q,cen); gc()
}

cat("3D nonconstant bin-pairs:", sum(Dgood), "/", nBP, "\n")

# ============================================================
# Randomized SVD for 3D
#
# Z is binpair x cell.
# Oversampled randomized range finder, with N_POWER power steps.
# ============================================================

L <- K + OVERSAMPLE

set.seed(RANDOM_SEED)
Omega <- matrix(rnorm(nC*L), nC, L)

cat("\n3D randomized SVD: initial range pass...\n")

Y <- matrix(0, nBP, L)

for (ib in seq_along(blocks)) {
    ii <- blocks[[ib]]
    z <- get_z3d(ii)
    Y[ii,] <- z %*% Omega

    if (ib%%10==0 || ib==length(blocks))
        cat("3D range",ib,"/",length(blocks),"\n")

    rm(z); gc()
}

Q <- qr.Q(qr(Y, LAPACK=TRUE))[,seq_len(L),drop=FALSE]
rm(Y,Omega); gc()

if (N_POWER > 0L) {

    for (qq in seq_len(N_POWER)) {

        cat("\n3D power iteration",qq,": Z'Q pass...\n")

        W <- matrix(0, nC, L)

        for (ib in seq_along(blocks)) {
            ii <- blocks[[ib]]
            z <- get_z3d(ii)
            W <- W + crossprod(z, Q[ii,,drop=FALSE])

            if (ib%%10==0 || ib==length(blocks))
                cat("3D power Z'Q",ib,"/",length(blocks),"\n")

            rm(z); gc()
        }

        cat("3D power iteration",qq,": ZW pass...\n")

        Y <- matrix(0, nBP, L)

        for (ib in seq_along(blocks)) {
            ii <- blocks[[ib]]
            z <- get_z3d(ii)
            Y[ii,] <- z %*% W

            if (ib%%10==0 || ib==length(blocks))
                cat("3D power ZW",ib,"/",length(blocks),"\n")

            rm(z); gc()
        }

        Q <- qr.Q(qr(Y, LAPACK=TRUE))[,seq_len(L),drop=FALSE]

        rm(W,Y); gc()
    }
}

# Small matrix Bsmall = Q'Z
cat("\n3D randomized SVD: Q'Z pass...\n")

Bsmall <- matrix(0, L, nC)

for (ib in seq_along(blocks)) {
    ii <- blocks[[ib]]
    z <- get_z3d(ii)
    Bsmall <- Bsmall + crossprod(Q[ii,,drop=FALSE], z)

    if (ib%%10==0 || ib==length(blocks))
        cat("3D Q'Z",ib,"/",length(blocks),"\n")

    rm(z); gc()
}

ss <- svd(Bsmall, nu=L, nv=L)

dD <- ss$d[seq_len(K)]
UD <- Q %*% ss$u[,seq_len(K),drop=FALSE]
VD <- ss$v[,seq_len(K),drop=FALSE]

# deterministic sign
for (k in seq_len(K)) {
    j <- which.max(abs(VD[,k]))
    if (VD[j,k] < 0) {
        VD[,k] <- -VD[,k]
        UD[,k] <- -UD[,k]
    }
}

normKD <- sqrt(sum(dD^2))
dD_scaled <- dD/normKD

BD <- sweep(UD, 2, dD_scaled, "*")

rownames(BD) <- as.character(bp)
colnames(BD) <- sprintf("K%02d",seq_len(K))
rownames(VD) <- final_cells
colnames(VD) <- colnames(BD)

total_fro2_D <- sum(pmax(Dn[Dgood]-1L,0L))
explD <- dD^2/total_fro2_D

# ============================================================
# Save
# ============================================================

outA <- list(
    modality="ATAC", K=K,
    feature_scores=SVD_A$B,
    cell_loadings=SVD_A$V,
    singular_values=SVD_A$d,
    explained_fraction=SVD_A$explained,
    retained_norm=SVD_A$normK
)

outH <- list(
    modality="H3K27ac", K=K,
    feature_scores=SVD_H$B,
    cell_loadings=SVD_H$V,
    singular_values=SVD_H$d,
    explained_fraction=SVD_H$explained,
    retained_norm=SVD_H$normK
)

outR <- list(
    modality="RNA", K=K,
    feature_scores=SVD_R$B,
    cell_loadings=SVD_R$V,
    singular_values=SVD_R$d,
    explained_fraction=SVD_R$explained,
    retained_norm=SVD_R$normK
)

outD <- list(
    modality="3D", K=K,
    feature_scores=BD,
    cell_loadings=VD,
    singular_values=dD,
    explained_fraction=explD,
    retained_norm=normKD,
    oversample=OVERSAMPLE,
    n_power=N_POWER,
    random_seed=RANDOM_SEED
)

stats <- list(
    cells=final_cells,
    enhancers=enh,
    genes=genes,
    binpairs=bp,
    ATAC_depth=ATAC_depth,
    H3K27ac_depth=H3_depth,
    ATAC_mean=SVD_A$stats$mu,
    ATAC_sd=SVD_A$stats$sd,
    H3_mean=SVD_H$stats$mu,
    H3_sd=SVD_H$stats$sd,
    RNA_mean=SVD_R$stats$mu,
    RNA_sd=SVD_R$stats$sd,
    D_mean=Dmu,
    D_sd=Dsd,
    D_nobs=Dn
)

saveRDS(outA, "GSE303006_M23_stage1_K20_ATAC.rds")
saveRDS(outH, "GSE303006_M23_stage1_K20_H3K27ac.rds")
saveRDS(outR, "GSE303006_M23_stage1_K20_RNA.rds")
saveRDS(outD, "GSE303006_M23_stage1_K20_3D.rds")
saveRDS(stats, "GSE303006_M23_stage1_K20_preprocessing_stats.rds")

# ---------- summaries ----------
cat("\n===== Stage 1 complete =====\n")
cat("ATAC feature scores:",paste(dim(outA$feature_scores),collapse=" x "),"\n")
cat("H3 feature scores  :",paste(dim(outH$feature_scores),collapse=" x "),"\n")
cat("RNA feature scores :",paste(dim(outR$feature_scores),collapse=" x "),"\n")
cat("3D feature scores  :",paste(dim(outD$feature_scores),collapse=" x "),"\n")

cat("\nATAC singular values:\n"); print(outA$singular_values)
cat("\nH3K27ac singular values:\n"); print(outH$singular_values)
cat("\nRNA singular values:\n"); print(outR$singular_values)
cat("\n3D singular values:\n"); print(outD$singular_values)

cat("\nStage-2 conceptual tensor:\n")
cat(nrow(ep)," x ",length(final_cells)," x ",K," x 4\n",sep="")

h5closeAll()
