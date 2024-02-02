# Compared weighted variance along *rows* of X
# where each element of X has its own positive weight
rowWeightedVarsMatrix <- function(X, W) {
  stopifnot(dim(X) == dim(W))

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(W)) W <- as.matrix(W)

  # scale weights to have a mean of 1
  Ws <- W / rowMeans2(W, useNames = FALSE)

  # weighted mean
  mu <- rowMeans2(X * Ws, useNames = FALSE)

  # weighted deviation
  A <- sqrt(Ws) * (X - mu)

  # sum of squares, divided by n-1
  rowSums2(A**2, useNames = FALSE) / (ncol(X) - 1)
}


#' @importFrom dplyr tibble bind_rows
getVarFromCounts <- function(countMatrix, lib.size, prior.count = .25) {
  stopifnot(ncol(countMatrix) == length(lib.size))

  # if lib.size is NA, set it to zero
  lib.size[is.na(lib.size)] <- 0

  countMatrix <- countMatrix + prior.count

  # pseudobulk
  count.gene <- rowSums2(countMatrix, useNames = FALSE)

  # normalize counts by library size
  # add pseudocount to counts here
  normCounts <- scale(countMatrix,
    scale = lib.size + 1,
    center = FALSE
  )

  # compute variance for each row
  # sigSq.mle <- rowVars(normCounts, useNames=FALSE)

  # weighted variance
  sigSq.mle <- rowWeightedVarsMatrix(normCounts, countMatrix)
  sigSq.mle[is.na(sigSq.mle)] <- 0

  # return values to compute variance later
  tibble(
    Gene = rownames(countMatrix),
    count.gene = count.gene,
    sigSq.mle = sigSq.mle,
    zeta = mean(lib.size^2),
    ncell = ncol(countMatrix)
  )
}


getVarForCellType <- function(sce, sample_id, cluster_id, geneList, CT, prior.count, verbose = TRUE) {
  cellType <- ID <- NULL

  if (!"counts" %in% assayNames(sce)) {
    stop("SCE does not contain assay: counts")
  }

  if (verbose) message("  Computing library sizes...")

  idx <- which(sce[[cluster_id]] == CT)

  if (is(counts(sce), "DelayedArray")) {
    lib.size <- colSums2(counts(sce), cols = idx, useNames =FALSE)
  } else {
    lib.size <- colSums2(counts(sce)[, idx], useNames = FALSE)
  }
  names(lib.size) = colnames(sce)[idx]

  # Scale prior count so that an observed count of 0,
  # gives zero variance across samples
  # Add small value to each cell, so that across n_i cells
  # it the augment sum to a mean of prior.count
  df_pc <- data.frame(
    ID = sce[[sample_id]][idx],
    cellType = sce[[cluster_id]][idx],
    prior.count = prior.count * lib.size / mean(lib.size)
  ) %>%
    group_by(cellType, ID) %>%
    summarize(n = length(ID), prior.count = mean(prior.count) / n, .groups = "drop_last")

  if (verbose) message("  Processing samples...")

  # genes tp keep
  keep <- geneList[[CT]]

  # get variance estimates for each ID and gene
  df <- lapply(unique(sce[[sample_id]]), function(ID) {

    idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == ID
    countMatrix <- counts(sce)[keep, idx, drop = FALSE]

    pc <- df_pc$prior.count[df_pc$ID == ID]

    if (ncol(countMatrix) > 0) {
      res <- getVarFromCounts(countMatrix,
        lib.size = lib.size[colnames(countMatrix)],
        prior.count = pc
      )
    } else {
      # if no cells are observed for this cell type and sample
      res <- tibble(
        Gene = keep,
        count.gene = 0,
        sigSq.mle = 0,
        zeta = 0,
        ncell = 0
      )
    }
    res$ID <- ID
    res
  })
  bind_rows(df)
}

#' @importFrom limma squeezeVar
#' @importFrom Matrix sparseMatrix
#' @importFrom dplyr mutate
getVarList <- function(sce, sample_id, cluster_id, geneList = NULL, shrink = TRUE, prior.count = 0.5, details = FALSE, verbose = TRUE) {
  Gene <- ID <- count.gene <- ncell <- zeta <- sigSq.mle <- vif <- NULL

  if (!sample_id %in% colnames(colData(sce))) {
    msg <- paste0("sample_id entry not found in colData(sce): ", sample_id)
    stop(msg)
  }
  if (!cluster_id %in% colnames(colData(sce))) {
    msg <- paste0("cluster_id entry not found in colData(sce): ", cluster_id)
    stop(msg)
  }

  # unique cluster ids
  clIDs <- as.character(unique(sce[[cluster_id]]))

  if( is.null(geneList) ){
    # create geneList with all genes
    geneList <- lapply(clIDs, function(x) rownames(sce))
    names(geneList) = clIDs
  }else{
    fnd = clIDs %in% names(geneList)
    if( ! all(fnd) ){
      msg <- paste(unique(sce[[cluster_id]])[!fnd], collapse=",")
      msg <- paste("assays not found in geneList:", msg)
      stop(msg)
    }

    # only include clusters with non-null gene sets
    exclude = names(geneList)[sapply(geneList, is.null)]
    clIDs = setdiff(clIDs, exclude)

    if( length(clIDs) == 0) stop("No data retained")
  }

  # Compute variance for each observation for each cell type
  var.list <- lapply(clIDs, function(CT) {
    if (verbose) message("Processing: ", CT)

    df <- getVarForCellType(sce, sample_id, cluster_id, geneList, CT, prior.count, verbose) %>%
      mutate(
        Gene = factor(Gene, geneList[[CT]]),
        ID = factor(ID)
      )

    if (shrink) {
      # shrink sample variances
      # use small offset to handle cases with zero variance
      # use pmax to ensure df > 1
      res <- squeezeVar(df$sigSq.mle + 1e-9, pmax(1, df$ncell - 1), robust = FALSE)
      df$sigSq.hat <- res$var.post
    }else{
      df$sigSq.hat <- df$sigSq.mle
    }

    # delta approximation of variance
    df <- df %>%
      mutate(count.gene = count.gene + 1e-4) %>%
      mutate(vif = (1 + ncell * sigSq.hat * zeta / count.gene)) %>%
      mutate(vhat = 1 / count.gene * vif) %>%
      mutate( assay = CT)

    # Vhat
    matVhat <- sparseMatrix(df$Gene, df$ID,
      x = df$vhat,
      dims = c(nlevels(df$Gene), nlevels(df$ID)),
      dimnames = list(levels(df$Gene), levels(df$ID))
    )
    matVhat <- as.matrix(matVhat)

    attr(matVhat, "details") <- df
    matVhat
  })
  names(var.list) <- clIDs

  var.list
}


# Get offset so that (max(x) + offset) / (min(x) + offset) is target_ratio
# return max(0, tau), so offset isn't negative
get_offset <- function(x, target_ratio) {
  # min and max
  rng <- range(x)

  # tau
  tau <- (rng[2] - target_ratio * rng[1]) / (target_ratio - 1)

  max(tau, 0)
}

#' Compute precision weights for pseudobulk
#'
#' Compute precision weights for pseudobulk using the delta method to approximate the variance of the log2 counts per million considering variation in the number of cells and gene expression variance across cells within each sample. By default, used number of cells; if specified use delta method.  Note that \code{processAssays()} uses number of cells as weights when no weights are specificed
#'
#' @param sce \code{SingleCellExperiment} of where \code{counts(sce)} stores the raw count data at the single cell level
#' @param sample_id character string specifying which variable to use as sample id
#' @param cluster_id character string specifying which variable to use as cluster id
#' @param geneList list of genes to be included for each cell type 
#' @param method select method to compute precision weights.  \code{'delta'} use the delta method based on normal approximation to a negative binomial model, slower but can increase power. \code{'ncells'} use the number of cells, this is faster; Subsequent arguments are ignored. Included for testing
#' @param shrink Defaults to \code{TRUE}. Use empirical Bayes variance shrinkage from \code{limma} to shrink estimates of expression variance across cells within each sample
#' @param prior.count Defaults to \code{0.5}. Count added to each observation at the pseudobulk level.  This is scaled but the number of cells before added to the cell level
# @param quantileOffset Defaults to \code{0.1}. When computing the precision from the variance, regularize the reciprocal by adding a small value to the denominator. For a gene with variances stored in the array \code{x}, add \code{quantile(x, quantileOffset)} before taking the reciprocal.
#' @param maxRatio When computing precision as the reciprocal of variance \code{1/(x+tau)} select tau to have a maximum ratio between the largest and smallest precision
#' @param h5adBlockSizes set the automatic block size block size (in bytes) for DelayedArray to read an H5AD file.  Larger values use more memory but are faster.
#' @param details include \code{data.frame} of cell-level statistics as \code{attr(., "details")}
#' @param verbose Show messages, defaults to TRUE
#'
#' @examples
#' library(muscat)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce,
#'   assay = "counts",
#'   sample_id = "sample_id",
#'   cluster_id = "cluster_id",
#'   verbose = FALSE
#' )
#'
#' # Gene expressed genes for each cell type
#' geneList = getExprGeneNames(pb)
#' 
#' # Create precision weights for pseudobulk
#' # By default, weights are set to cell count,
#' # which is the default in processAssays()
#' # even when no weights are specified
#' weightsList <- pbWeights(example_sce,
#'   sample_id = "sample_id",
#'   cluster_id = "cluster_id",
#'   geneList = geneList
#' )
#'
#' # voom-style normalization using initial weights
#' res.proc <- processAssays(pb, ~group_id, weightsList = weightsList)
#
#' @importFrom stats quantile
#' @importFrom DelayedArray getAutoBlockSize setAutoBlockSize
#' @export
pbWeights <- function(sce, sample_id, cluster_id, geneList = NULL, method = c("delta", "ncells"), shrink = TRUE, prior.count = 0.5, maxRatio = 20, h5adBlockSizes = 1e9, details = FALSE, verbose = TRUE) {
  method <- match.arg(method)

  # check for NA values
  if (any(is.na(sce[[sample_id]]))) {
    stop("NA values are not allowed in sample_id column")
  }
  if (any(is.na(sce[[cluster_id]]))) {
    stop("NA values are not allowed in cluster_id column")
  }

  if (method == "ncells") {
    W.list <- .pbWeights_ncells(sce, sample_id, cluster_id)
  } else {
    # delta approximation

    # update block size for reading h5ad file from disk
    tmp <- getAutoBlockSize()
    suppressMessages(setAutoBlockSize(h5adBlockSizes))
    on.exit(suppressMessages(setAutoBlockSize(tmp)))

    # compute variances
    var.lst <- getVarList(sce, sample_id, cluster_id, geneList, shrink, prior.count, details = details, verbose = verbose)

    # for each cell type
    W.list <- lapply( names(var.lst), function(id) {
      v.mat = var.lst[[id]]
      # regularize reciprocal with offset
      # get offset tau as the max of x
      #  to give a maximum ratio of maxRatio
      # for each cell type
      ret <- t(apply(v.mat, 1, function(x) {
        tau <- get_offset(x, maxRatio)
        1 / (x + tau)
      }))

      if( details ) attr(ret, "details") <- attr(var.lst[[id]], "details") 
      ret
    })
    names(W.list) <- names(var.lst)
  }

  W.list
}


.pbWeights_ncells <- function(sce, sample_id, cluster_id) {
  if (!sample_id %in% colnames(colData(sce))) {
    txt <- paste("sample_id not found:", sample_id)
    stop(txt)
  }
  if (!cluster_id %in% colnames(colData(sce))) {
    txt <- paste("cluster_id not found:", cluster_id)
    stop(txt)
  }

  # number of cells
  df_cc <- table(sce[[sample_id]], sce[[cluster_id]])

  W.list <- lapply(unique(sce[[cluster_id]]), function(k) {
    W <- matrix(df_cc[, k], ncol = nrow(df_cc), nrow = nrow(sce), byrow = TRUE)
    colnames(W) <- rownames(df_cc)
    rownames(W) <- rownames(sce)
    W
  })
  names(W.list) <- unique(sce[[cluster_id]])
  W.list
}
