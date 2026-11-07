# ============================================================
# revision_lib.R
# GSE303006 / GENCODE M23
#
# Shared helper functions for the revision analyses (scripts 19-27).
# This file only defines functions; it does not compute anything.
#
# Contents
#   rev_find_repo()            locate the repository (for sourcing)
#   rev_start_log()/rev_end_log()
#                              save the console output of a run
#   rev_kw_epsilon2()          Kruskal-Wallis epsilon^2 (same as 10)
#   rev_rank_biserial()        Wilcoxon rank-biserial effect
#   rev_colorbar()             vertical colour bar for base-graphics
#                              heatmaps
#   rev_stage2_hosvd()         implicit Stage-2 HOSVD for an arbitrary
#                              subset of modalities (same algebra as 08)
#   rev_ep_factor_column()     one column of U_EP, computed blockwise
#   rev_contrast_groups()      resolve the two pre-specified cell types
#   rev_select_core()          pre-specified core-selection rule
#   rev_gene_scores()          degree-corrected gene score (same as 11)
#   rev_run_gsea()             preranked GO-BP GSEA (same settings as 12)
#   rev_compare_rankings()     agreement of two gene rankings
#   rev_compare_go()           agreement of two GSEA results
# ============================================================

suppressPackageStartupMessages(library(data.table))


# ------------------------------------------------------------
# Repository location
#
# Priority: environment variable TDCHARM_REPO, then the variable
# `repo` set by run_analysis.R, then the directory of the script
# being sourced, then the working directory.
# ------------------------------------------------------------

rev_find_repo <- function() {

    env <- Sys.getenv("TDCHARM_REPO", unset = "")
    if (nzchar(env) && dir.exists(file.path(env, "R"))) {
        return(normalizePath(env))
    }

    if (exists("repo", envir = globalenv(), inherits = FALSE)) {
        r <- get("repo", envir = globalenv())
        if (is.character(r) && dir.exists(file.path(r, "R"))) {
            return(normalizePath(r))
        }
    }

    for (i in rev(seq_len(sys.nframe()))) {
        of <- sys.frame(i)$ofile
        if (!is.null(of)) {
            d <- dirname(dirname(normalizePath(of)))
            if (dir.exists(file.path(d, "R"))) return(d)
        }
    }

    args <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args, value = TRUE)
    if (length(fa)) {
        d <- dirname(dirname(normalizePath(sub("^--file=", "", fa[1]))))
        if (dir.exists(file.path(d, "R"))) return(d)
    }

    if (dir.exists("R")) return(normalizePath("."))

    stop("Cannot locate the TDCHARM repository; set TDCHARM_REPO.")
}


# ------------------------------------------------------------
# Logging: everything printed to the console is also written to
# the log file. The log is stored next to the outputs.
# ------------------------------------------------------------

rev_start_log <- function(log_file, script_name) {
    con <- file(log_file, open = "wt")
    sink(con, split = TRUE)
    cat(script_name, "\n", sep = "")
    cat("Run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
    cat("Working directory: ", getwd(), "\n", sep = "")
    cat("R: ", R.version.string, "\n\n", sep = "")
    invisible(con)
}

rev_end_log <- function(con) {
    cat("\nSession info:\n")
    print(sessionInfo())
    sink()
    close(con)
    invisible(NULL)
}


# ------------------------------------------------------------
# Effect sizes
# ------------------------------------------------------------

# Identical definition to kw_epsilon2() in script 10.
rev_kw_epsilon2 <- function(x, group) {
    ok <- is.finite(x) & !is.na(group)
    x <- x[ok]
    group <- factor(group[ok])
    if (length(x) < 2L || nlevels(group) < 2L) {
        return(list(statistic = NA_real_, df = NA_real_,
                    p = NA_real_, epsilon2 = NA_real_, n = length(x)))
    }
    kt <- kruskal.test(x ~ group)
    H <- as.numeric(kt$statistic)
    klev <- nlevels(group)
    n <- length(x)
    eps2 <- (H - klev + 1) / (n - klev)
    eps2 <- max(0, min(1, eps2))
    list(statistic = H, df = as.numeric(kt$parameter),
         p = kt$p.value, epsilon2 = eps2, n = n)
}

# Rank-biserial correlation of x versus y (positive: x > y).
rev_rank_biserial <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (!length(x) || !length(y)) return(NA_real_)
    r <- rank(c(x, y))
    U <- sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2
    2 * U / (length(x) * length(y)) - 1
}


# ------------------------------------------------------------
# Vertical colour bar drawn in the right margin of the current
# plot (base graphics).
# ------------------------------------------------------------

rev_colorbar <- function(cols, zlim, title = "", at = NULL,
                         x_from = 1.02, x_to = 1.05) {
    usr <- par("usr")
    xr <- usr[2] - usr[1]
    x0 <- usr[2] + (x_from - 1) * xr
    x1 <- usr[2] + (x_to - 1) * xr
    n <- length(cols)
    y <- seq(usr[3], usr[4], length.out = n + 1)
    op <- par(xpd = NA)
    on.exit(par(op))
    rect(x0, y[-(n + 1)], x1, y[-1], col = cols, border = NA)
    rect(x0, usr[3], x1, usr[4], border = "black")
    if (is.null(at)) at <- pretty(zlim)
    at <- at[at >= zlim[1] & at <= zlim[2]]
    ya <- usr[3] + (at - zlim[1]) / diff(zlim) * (usr[4] - usr[3])
    segments(x1, ya, x1 + 0.006 * xr, ya)
    text(x1 + 0.010 * xr, ya, labels = format(at), adj = c(0, 0.5),
         cex = 0.75)
    text((x0 + x1) / 2, usr[4] + 0.03 * (usr[4] - usr[3]), title,
         adj = c(0.5, 0), cex = 0.8)
    invisible(NULL)
}


# ------------------------------------------------------------
# Implicit Stage-2 HOSVD
#
# Same algebra as 08_stage2_implicit_HOSVD.R, generalized to any
# subset of the four modalities (for the modality-ablation and
# rotation analyses).
#
#   S          named list of Stage-1 objects (feature_scores,
#              cell_loadings), one per modality, in the order used
#              for the modality mode
#   map        Stage-2 mapping table from 07 (idx_*, w_*)
#   R_EP, R_CELL
#              requested ranks (capped at the numerical rank)
#
#   X[p,c,k,m] = w_m[p] B_m[idx_m(p),k] V_m[c,k]
#
#   E-P factor   : eigenvectors of H_A = A^T A (A = pair x (M K))
#   cell factor  : left singular vectors of [ ||a_mk|| v_mk ]
#   K factor     : identity (K-mode Gram is diagonal)
#   modality     : eigenvectors of G_M[m,n] =
#                    sum_k <a_mk,a_nk><v_mk,v_nk>
#   core         : CORE[a,b,k,r] =
#                    sum_m PAIR_PROJ[a,(m,k)] CELL_PROJ[b,(m,k)] U_MODAL[m,r]
#
# The full E-P factor matrix is not returned (it is P x R_EP);
# rev_ep_factor_column() rebuilds any column blockwise.
# ------------------------------------------------------------

REV_IDX_COL <- c(ATAC = "idx_ATAC", H3K27ac = "idx_H3K27ac",
                 RNA = "idx_RNA", `3D` = "idx_3D")
REV_W_COL <- c(ATAC = "w_enhancer", H3K27ac = "w_enhancer",
               RNA = "w_gene", `3D` = "w_binpair")

rev_build_A_block <- function(S, map, ii) {
    mp <- map[ii]
    do.call(cbind, lapply(names(S), function(m) {
        B <- S[[m]]$feature_scores
        B[mp[[REV_IDX_COL[[m]]]], , drop = FALSE] * mp[[REV_W_COL[[m]]]]
    }))
}

rev_stage2_hosvd <- function(S, map, R_EP = 20L, R_CELL = 20L,
                             pair_block = 10000L, verbose = TRUE) {

    mods <- names(S)
    stopifnot(all(mods %in% names(REV_IDX_COL)))
    M <- length(mods)
    K <- ncol(S[[1]]$feature_scores)
    cells <- rownames(S[[1]]$cell_loadings)
    for (m in mods) {
        stopifnot(ncol(S[[m]]$feature_scores) == K,
                  identical(rownames(S[[m]]$cell_loadings), cells))
    }
    P <- nrow(map)
    L <- M * K
    offset <- lapply(seq_len(M), function(i) ((i - 1L) * K + 1L):(i * K))
    names(offset) <- mods

    blocks <- split(seq_len(P), ceiling(seq_len(P) / pair_block))

    # Pair-side Gram matrix
    H_A <- matrix(0, L, L)
    for (ib in seq_along(blocks)) {
        Ab <- rev_build_A_block(S, map, blocks[[ib]])
        H_A <- H_A + crossprod(Ab)
    }
    H_A <- (H_A + t(H_A)) / 2

    modal_norm2 <- vapply(offset, function(jj) sum(diag(H_A)[jj]), numeric(1))
    tensor_fro2 <- sum(modal_norm2)

    # E-P mode
    eeP <- eigen(H_A, symmetric = TRUE)
    lambdaP_all <- pmax(eeP$values, 0)
    rankP <- sum(lambdaP_all > max(lambdaP_all) * 1e-12)
    R_EP <- min(R_EP, rankP)
    WP <- eeP$vectors[, seq_len(R_EP), drop = FALSE]
    sigmaP <- sqrt(lambdaP_all[seq_len(R_EP)])

    # Cell mode
    V <- do.call(cbind, lapply(mods, function(m) S[[m]]$cell_loadings))
    normA <- sqrt(pmax(diag(H_A), 0))
    Ccat <- sweep(V, 2, normA, "*")
    svC <- svd(Ccat, nu = min(nrow(Ccat), ncol(Ccat)), nv = 0)
    R_CELL <- min(R_CELL, sum(svC$d > max(svC$d) * 1e-12))
    U_CELL <- svC$u[, seq_len(R_CELL), drop = FALSE]
    rownames(U_CELL) <- cells
    colnames(U_CELL) <- sprintf("CELL%02d", seq_len(R_CELL))

    # Modality mode
    G_M <- matrix(0, M, M)
    for (i in seq_len(M)) for (j in seq_len(M)) {
        pc <- diag(H_A[offset[[i]], offset[[j]], drop = FALSE])
        cc <- diag(crossprod(S[[mods[i]]]$cell_loadings,
                             S[[mods[j]]]$cell_loadings))
        G_M[i, j] <- sum(pc * cc)
    }
    G_M <- (G_M + t(G_M)) / 2
    eeM <- eigen(G_M, symmetric = TRUE)
    U_MODAL <- eeM$vectors
    rownames(U_MODAL) <- mods
    colnames(U_MODAL) <- sprintf("MOD%02d", seq_len(M))

    # Core
    PAIR_PROJ <- diag(sigmaP, nrow = R_EP) %*% t(WP)        # R_EP x L
    core_pre <- array(0, dim = c(R_EP, R_CELL, K, M))
    for (i in seq_len(M)) {
        Cp <- crossprod(U_CELL, S[[mods[i]]]$cell_loadings)  # R_CELL x K
        jj <- offset[[i]]
        for (k in seq_len(K)) {
            core_pre[, , k, i] <- tcrossprod(PAIR_PROJ[, jj[k]], Cp[, k])
        }
    }
    CORE <- array(0, dim = c(R_EP, R_CELL, K, M),
                  dimnames = list(sprintf("EP%02d", seq_len(R_EP)),
                                  colnames(U_CELL),
                                  sprintf("K%02d", seq_len(K)),
                                  colnames(U_MODAL)))
    for (r in seq_len(M)) for (i in seq_len(M)) {
        CORE[, , , r] <- CORE[, , , r] + core_pre[, , , i] * U_MODAL[i, r]
    }

    core_fro2 <- sum(CORE^2)
    if (verbose) {
        cat(sprintf("Stage 2 [%s]: P=%d, cells=%d, K=%d, M=%d; ranks %d x %d x %d x %d; core energy fraction = %.4f\n",
                    paste(mods, collapse = "+"), P, length(cells), K, M,
                    R_EP, R_CELL, K, M, core_fro2 / tensor_fro2))
    }

    list(modalities = mods, K = K, cells = cells, offset = offset,
         H_A = H_A, WP = WP, sigmaP = sigmaP, U_CELL = U_CELL,
         U_MODAL = U_MODAL, CORE = CORE, PAIR_PROJ = PAIR_PROJ,
         modal_norm2 = modal_norm2, tensor_fro2 = tensor_fro2,
         core_fro2 = core_fro2, lambdaP_all = lambdaP_all,
         sigmaC_all = svC$d, G_M = G_M)
}

# Column a of U_EP:  U_EP = A W Sigma^{-1}
rev_ep_factor_column <- function(fit, S, map, a, pair_block = 10000L) {
    P <- nrow(map)
    blocks <- split(seq_len(P), ceiling(seq_len(P) / pair_block))
    u <- numeric(P)
    for (ii in blocks) {
        u[ii] <- rev_build_A_block(S, map, ii) %*% fit$WP[, a] / fit$sigmaP[a]
    }
    u
}


# ------------------------------------------------------------
# Pre-specified contrast cell types
#
# Defaults correspond to the contrast reported for core #12
# (Inh Ndnf/Lamp5 versus Ex L3/4 IT). They are resolved against
# the cell-type labels of script 10 by regular expression and can
# be overridden with environment variables containing exact
# labels separated by ';':
#   TDCHARM_GROUP_A, TDCHARM_GROUP_B
# ------------------------------------------------------------

REV_GROUP_A_REGEX <- "Ndnf|Lamp5"
REV_GROUP_B_REGEX <- "L3.?4.*IT"

rev_contrast_groups <- function(cell_types) {
    labs <- sort(unique(na.omit(cell_types)))
    pick <- function(env, regex) {
        e <- Sys.getenv(env, unset = "")
        if (nzchar(e)) {
            x <- trimws(strsplit(e, ";", fixed = TRUE)[[1]])
            bad <- setdiff(x, labs)
            if (length(bad)) stop(env, ": unknown label(s): ", paste(bad, collapse = ", "))
            return(x)
        }
        x <- grep(regex, labs, value = TRUE)
        if (!length(x)) {
            stop("No cell-type label matches /", regex, "/. Available labels:\n",
                 paste(labs, collapse = "\n"),
                 "\nSet ", env, " to the exact label(s).")
        }
        x
    }
    A <- pick("TDCHARM_GROUP_A", REV_GROUP_A_REGEX)
    B <- pick("TDCHARM_GROUP_B", REV_GROUP_B_REGEX)
    if (length(intersect(A, B))) stop("Contrast groups overlap.")
    list(A = A, B = B)
}


# ------------------------------------------------------------
# Pre-specified core-selection rule
#
# Mirrors how core #12 / CELL03 is characterized in the manuscript
# (cell types enriched among the 50 most positive and 50 most
# negative cells), without using any decomposition other than the
# one being analysed.
#
#   1. For every cell component b, take the TOP_N cells at each end.
#      score_b = max over the two orientations of
#                frac(A among one end) + frac(B among the other end)
#      Ties are broken by |rank-biserial(A vs B)| over all cells.
#   2. b* = argmax score_b.
#   3. (a*, k*, r*) = argmax_{a,k,r} |CORE[a, b*, k, r]|.
#   4. Orientation: the modality column r* is oriented to a
#      positive sum and the cell column b* so that group A has the
#      larger mean loading. The oriented core value is
#         G* = CORE[a*,b*,k*,r*] * s_r * s_b
#      and the E-P score s_p = G* U_EP[p,a*] is invariant to the
#      arbitrary sign of U_EP[,a*].
# ------------------------------------------------------------

rev_select_core <- function(U_CELL, CORE, U_MODAL, cell_type, groups,
                            top_n = 50L) {

    stopifnot(length(cell_type) == nrow(U_CELL))
    inA <- cell_type %in% groups$A
    inB <- cell_type %in% groups$B

    comp <- rbindlist(lapply(seq_len(ncol(U_CELL)), function(b) {
        x <- U_CELL[, b]
        o <- order(x, decreasing = TRUE)
        hi <- o[seq_len(top_n)]
        lo <- rev(o)[seq_len(top_n)]
        s1 <- mean(inA[hi]) + mean(inB[lo])   # A high, B low
        s2 <- mean(inA[lo]) + mean(inB[hi])   # A low, B high
        data.table(b = b, CELL_component = colnames(U_CELL)[b],
                   score = max(s1, s2),
                   fracA_end = if (s1 >= s2) mean(inA[hi]) else mean(inA[lo]),
                   fracB_end = if (s1 >= s2) mean(inB[lo]) else mean(inB[hi]),
                   rank_biserial_A_vs_B = rev_rank_biserial(x[inA], x[inB]),
                   nA = sum(inA), nB = sum(inB))
    }))
    setorder(comp, -score, -rank_biserial_A_vs_B)
    comp[, abs_rb := abs(rank_biserial_A_vs_B)]
    setorder(comp, -score, -abs_rb)
    b <- comp$b[1]

    sub <- CORE[, b, , , drop = FALSE]
    oo <- order(abs(as.numeric(sub)), decreasing = TRUE)
    idx_all <- arrayInd(oo[seq_len(min(5L, length(oo)))], dim(sub))
    candidates <- data.table(
        candidate = seq_len(nrow(idx_all)), EP = idx_all[, 1], CELL = b,
        K = idx_all[, 3], MOD = idx_all[, 4],
        abs_core = abs(as.numeric(sub)[oo[seq_len(nrow(idx_all))]]))
    candidates[, ratio_to_selected := abs_core / abs_core[1]]
    idx <- idx_all[1, , drop = FALSE]
    a <- idx[1]; k <- idx[3]; r <- idx[4]

    s_b <- if (mean(U_CELL[inA, b]) >= mean(U_CELL[inB, b])) 1 else -1
    s_r <- if (sum(U_MODAL[, r]) >= 0) 1 else -1
    G_raw <- CORE[a, b, k, r]
    G_oriented <- G_raw * s_b * s_r

    core_rank <- sum(abs(as.numeric(CORE)) > abs(G_raw)) + 1L

    list(a = a, b = b, k = k, r = r, G_raw = G_raw,
         G_oriented = G_oriented, sign_cell = s_b, sign_modal = s_r,
         core_rank = core_rank, component_scores = comp,
         core_candidates = candidates,
         cell_loading_oriented = U_CELL[, b] * s_b)
}


# ------------------------------------------------------------
# Degree-corrected gene score (identical to script 11)
#
#   s_p = G * U_EP[p,a]
#   S_g = (1/sqrt(d_g)) sum_{p: g(p)=g} s_p
# ------------------------------------------------------------

rev_gene_scores <- function(map, ep_score) {
    stopifnot(length(ep_score) == nrow(map))
    dt <- data.table(gene_name = map$gene_name, d_g = map$d_g, s = ep_score)
    out <- dt[, .(gene_score = sum(s) / sqrt(d_g[1]), n_EP = .N),
              by = gene_name]
    setorder(out, -gene_score)
    out[]
}


# ------------------------------------------------------------
# Preranked GO-BP GSEA, same settings as script 12
# (clusterProfiler::gseGO, SYMBOL keys, sizes 10-500, BH, seed).
# ------------------------------------------------------------

rev_run_gsea <- function(gene_dt, seed = 12345L, min_gs = 10L,
                         max_gs = 500L) {
    for (p in c("clusterProfiler", "org.Mm.eg.db")) {
        if (!requireNamespace(p, quietly = TRUE)) stop("Package required: ", p)
    }
    gl <- gene_dt$gene_score
    names(gl) <- gene_dt$gene_name
    gl <- gl[is.finite(gl)]
    gl <- gl[!duplicated(names(gl))]
    gl <- sort(gl, decreasing = TRUE)
    # gseGO -> fgsea runs through BiocParallel, whose default back end
    # forks one worker per core of the NODE (not of the job). On a
    # 56-core node every worker copies the R session, and the job is
    # killed for exceeding its memory limit. Run GSEA serially.
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
        BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
    }
    options(mc.cores = 1L)
    set.seed(seed)
    res <- clusterProfiler::gseGO(
        geneList = gl, ont = "BP", OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        keyType = "SYMBOL", exponent = 1, minGSSize = min_gs,
        maxGSSize = max_gs, pvalueCutoff = 1, pAdjustMethod = "BH",
        verbose = FALSE)
    rev_gsea_table(as.data.frame(res))
}


# GSEA result -> data.table ordered by adjusted P, then |NES|.
# (data.table::setorder() accepts column names only, not expressions.)
rev_gsea_table <- function(df) {
    out <- as.data.table(df)
    if (!nrow(out)) return(out)
    out <- out[order(p.adjust, -abs(NES))]
    out[]
}

# ------------------------------------------------------------
# Agreement between two gene rankings (same gene universe).
#   Spearman correlation of signed scores, and overlap of the
#   top-N genes on the positive and negative sides with a
#   hypergeometric P value (universe = common genes).
# ------------------------------------------------------------

rev_compare_rankings <- function(x, y, label_x, label_y,
                                 top_n = c(100L, 500L)) {
    m <- merge(x[, .(gene_name, sx = gene_score)],
               y[, .(gene_name, sy = gene_score)], by = "gene_name")
    U <- nrow(m)
    rho <- suppressWarnings(cor(m$sx, m$sy, method = "spearman"))
    top_n <- unique(pmin(top_n, floor(U / 2)))
    rbindlist(lapply(top_n, function(n) {
        rbindlist(lapply(c("positive", "negative"), function(side) {
            sg <- if (side == "positive") 1 else -1
            gx <- m$gene_name[order(-sg * m$sx)][seq_len(n)]
            gy <- m$gene_name[order(-sg * m$sy)][seq_len(n)]
            ov <- length(intersect(gx, gy))
            data.table(run_x = label_x, run_y = label_y, n_genes = U,
                       spearman = rho, side = side, top_n = n,
                       overlap = ov, expected = n * n / U,
                       fold = ov / (n * n / U),
                       hypergeom_p = phyper(ov - 1, n, U - n, n,
                                            lower.tail = FALSE))
        }))
    }))
}


# ------------------------------------------------------------
# Agreement between two GSEA results.
#   For each side (NES>0, NES<0): number of terms at FDR <= cut in
#   each run, their overlap and Jaccard index; Spearman correlation
#   of NES over terms tested in both runs.
# ------------------------------------------------------------

rev_compare_go <- function(gx, gy, label_x, label_y, fdr = 0.05) {
    if (!nrow(gx) || !nrow(gy)) return(data.table())
    m <- merge(gx[, .(ID, NESx = NES)], gy[, .(ID, NESy = NES)], by = "ID")
    rho <- suppressWarnings(cor(m$NESx, m$NESy, method = "spearman"))
    rbindlist(lapply(c("positive", "negative"), function(side) {
        sx <- if (side == "positive") gx[p.adjust <= fdr & NES > 0, ID]
              else gx[p.adjust <= fdr & NES < 0, ID]
        sy <- if (side == "positive") gy[p.adjust <= fdr & NES > 0, ID]
              else gy[p.adjust <= fdr & NES < 0, ID]
        ov <- length(intersect(sx, sy))
        un <- length(union(sx, sy))
        data.table(run_x = label_x, run_y = label_y, side = side,
                   n_sig_x = length(sx), n_sig_y = length(sy),
                   overlap = ov, jaccard = if (un) ov / un else NA_real_,
                   NES_spearman_common_terms = rho,
                   n_common_terms = nrow(m))
    }))
}

# Top terms per side, without omission (for supplementary tables).
rev_go_top_terms <- function(g, label, n = 10L) {
    if (!nrow(g)) return(data.table())
    rbindlist(lapply(c("positive", "negative"), function(side) {
        z <- if (side == "positive") g[NES > 0] else g[NES < 0]
        z <- z[order(p.adjust, -abs(NES))]
        z <- z[seq_len(min(n, .N)),
               .(ID, Description, setSize, NES, pvalue, p.adjust)]
        z[, `:=`(run = label, side = side, rank = seq_len(.N))]
        setcolorder(z, c("run", "side", "rank"))
        z
    }))
}
