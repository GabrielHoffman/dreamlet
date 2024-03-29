# Gabriel Hoffman
# Nov 11, 2021



# min.count = 5
# normalize.method = "TMM"
# method = c("random", "fixed", "none"),



#' Differential expression between pair of assays
#'
#' Perform differential expression between a pair of assays using linear (mixed) models
#'
#' @param pb pseudobulk data as SingleCellExperiment object
#' @param assays array of two entries specifying assays (i.e. cell clusters) to compare, or a list of two sets of assays.
#' @param method account for repeated measures from donors using a "random" effect, a "fixed" effect, or "none"
#' @param formula covariates to include in the analysis.
#' @param collapse if TRUE (default), combine all cell clusters within the test set, and separately the baseline set. If FALSE, estimate coefficient for each cell cluster and then identify differential expression using linear contrasts with \code{variancePartition::makeContrastsDream()}
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count minimum number of reads for a gene to be consider expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param robust logical, use eBayes method that is robust to outlier genes
#' @param quiet show messages
#' @param contrasts cell type is encoded in variable \code{cellCluster} with levels \code{test} and \code{baseline}. \code{contrasts} specifies contrasts passed to \code{variancePartition::makeContrastsDream()}.  Note, advanced users only.
#' @param BPPARAM parameters for parallel evaluation
#' @param errorsAsWarnings if \code{TRUE}, convert error to a warning and return \code{NULL}
#' @param ... other arguments passed to \code{dream}
#'
#' @details
#' Analyze pseudobulk data to identify differential gene expression between two cell clusters or sets of clusters while modeling the cross-donor expression variation and other aspects of the study design.
#'
#' \code{dreamletCompareClusters()} is useful for finding genes that are differentially expressed betweeen cell clusters and estimating their fold change.  However, the p-values and number of differentially expressed genes are problematic for two reasons, so users must be careful not to overinterpret them:
#'
#' \enumerate{
#' \item{Cell clusters are typically identified with the same gene expression data used for this differential expression analysis between clusters. The same data is used both for discovery and testing, and this means that the p-values from the differential expression analysis will not be uniform under the null.  This will produce a lot of findings with small p-values even in the absence of true biological differences.}
#' \item{The \code{dreamlet} package is designed for large datasets with many subjects.  The sample sizes from cohort studies are an order of magnitude larger than typical single cell studies.  This means that these analyses have huge power to detect even subtle difference in expression between cell clusters.  While cluster-specific marker genes are often discovered from an handful of samples, the \code{dreamlet} package is applicable to 100s or 1000s of subjects.}
#' }
#'
#' \code{method} indicates the regression method used to test differential expression between sets of cell clusters.  Since the same biosample will usually be represented in both sets of cell clusters, \code{method} determines how the paired design is modeled.   For \code{method = "mixed"}, the sample is modeled as a random effect: \code{~ (1|Sample) + ...}. For \code{method = "fixed"}, the sample is modeled as a fixed effect: \code{~ Sample + ...}. For \code{method = "none"}, the pairing is ignored.
#'
#' When \code{collapse=TRUE} (default) combine all cell clusters within the test set, and separately the baseline set, and estimate a coefficient indicating the differential expression between sets for a given gene.  If \code{collapse=FALSE}, estimate a coefficient for each cell type and then identify differential expression using linear contrasts with \code{variancePartition::makeContrastsDream()}.
#'
#' @return Object of class \code{dreamletResult} storing results for each comparison
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce,
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#'
#' # Evaluate the specificity of each gene for each cluster
#' df_cts <- cellTypeSpecificity(pb)
#'
#' # compare first two assays (i.e. cell types)
#' ct.pairs <- c("B cells", "CD14+ Monocytes")
#'
#' # run comparison
#' # use method = 'fixed' here since it is faster
#' fit <- dreamletCompareClusters(pb, ct.pairs, method = "fixed")
#'
#' # Extract top 10 differentially expressed genes
#' # The coefficient 'compare' is the value logFC between test and baseline:
#' # compare = cellClustertest - cellClusterbaseline
#' res <- topTable(fit, coef = "compare", number = 10)
#'
#' # genes with highest logFC are most highly expressed in
#' # B cells compared to CD14+ Monocytes
#' head(res)
#'
#' dreamlet::plotHeatmap(df_cts, genes = rownames(res)[1:5])
#'
#' # compare B cells versus the rest of the cell types
#' # 'rest' is a keyword indicating all other assays
#' fit <- dreamletCompareClusters(pb, c("B cells", "rest"), method = "fixed")
#'
#' res <- topTable(fit, coef = "compare", number = 10)
#'
#' # genes with highest logFC are most highly expressed in
#' # B cells compared to all others
#' head(res)
#'
#' # Get genes upregulated in B cells
#' idx <- with(res, which(logFC > 0))[1:5]
#' dreamlet::plotHeatmap(df_cts, genes = rownames(res)[idx])
#'
#' lst <- list(
#'   test = c("CD14+ Monocytes", "FCGR3A+ Monocytes"),
#'   baseline = c("CD4 T cells", "CD8 T cells")
#' )
#'
#' # compare 2 monocyte clusters to two T cell clusters
#' fit <- dreamletCompareClusters(pb, lst, method = "fixed")
#'
#' res <- topTable(fit, coef = "compare", number = 10)
#'
#' # genes with highest logFC are most highly expressed in
#' # monocytes compared to T cells
#' head(res)
#'
#' # Get genes upregulated in monocytes
#' idx <- with(res, which(logFC > 0))[1:5]
#' dreamlet::plotHeatmap(df_cts, genes = rownames(res)[idx])
#'
#' @importFrom variancePartition dream eBayes topTable makeContrastsDream
#' @export
dreamletCompareClusters <- function(pb, assays, method = c("fixed", "random", "none"), formula = ~0, collapse = TRUE, min.cells = 10, min.count = 10, min.samples = 4, isCounts = TRUE, normalize.method = "TMM", robust = FALSE, quiet = FALSE, contrasts = c(compare = paste("cellClustertest - cellClusterbaseline")), BPPARAM = SerialParam(), errorsAsWarnings = FALSE, ...) {
  method <- match.arg(method)
  formula <- as.formula(formula)

  if (!is(pb, "SingleCellExperiment")) {
    stop("pb must be a SingleCellExperiment")
  }

  if (!all(all.vars(formula) %in% colnames(get_metadata_aggr_means(pb)))) {
    i <- which(!all.vars(formula) %in% colnames(get_metadata_aggr_means(pb)))
    cmd <- paste0("Variables in formula not found: ", paste(all.vars(formula)[i], collapse = ", "))
    stop(cmd)
  }

  if (is.vector(assays) & !is.list(assays)) {
    # check that assays has two entries
    if (length(assays) != 2) {
      txt <- "assays must have 2 entries"
      if (errorsAsWarnings) {
        warning(txt)
        return(NULL)
      } else {
        stop(txt)
      }
    }

    if (assays[1] == "rest") {
      txt <- "assay 'rest' is only valid as the second entry"
      if (errorsAsWarnings) {
        warning(txt)
        return(NULL)
      } else {
        stop(txt)
      }
    }

    # convert array to list with two entries
    if ("rest" %in% assays) {
      j <- which(assayNames(pb) == assays[1])
      assay.lst <- list(test = assays[1], baseline = assayNames(pb)[-j])
    } else {
      assay.lst <- list(test = assays[1], baseline = assays[2])
    }
  } else if (is.list(assays)) {
    assay.lst <- assays
  } else {
    txt <- "assays argument must be vector or list"
    if (errorsAsWarnings) {
      warning(txt)
      return(NULL)
    } else {
      stop(txt)
    }
  }

  # check that assay.lst has two entries
  if (length(assay.lst) != 2) {
    txt <- "assay.lst must have 2 entries"
    if (errorsAsWarnings) {
      warning(txt)
      return(NULL)
    } else {
      stop(txt)
    }
  }

  # extract metadata shared across assays
  data_constant <- as.data.frame(colData(pb))

  # remove samples with missing covariate data
  idx <- lapply(all.vars(formula), function(v) {
    which(is.na(data_constant[[v]]))
  })
  idx <- unique(unlist(idx))

  if (length(idx) > 1) {
    data_constant <- droplevels(data_constant[-idx, , drop = FALSE])
  }

  # convert entries to strings
  assay.lst <- lapply(assay.lst, as.character)

  # check that all are valid assays
  if (any(!unlist(assay.lst) %in% assayNames(pb))) {
    i <- unlist(assay.lst) %in% assayNames(pb)

    txt <- paste(unlist(assay.lst)[!i], collapse = ", ")
    txt <- paste("Specified assays are not valid: ", txt)

    if (errorsAsWarnings) {
      warning(txt)
      return(NULL)
    } else {
      stop(txt)
    }
  }

  # check that there are no overlapping assays
  shared <- intersect(assay.lst[[1]], assay.lst[[2]])
  if (length(shared) > 0) {
    txt <- paste(shared, collapse = ", ")
    txt <- paste("Specified assays shared between two groups: ", txt)

    if (errorsAsWarnings) {
      warning(txt)
      return(NULL)
    } else {
      stop(txt)
    }
  }

  if (collapse) {
    # for test and baseline
    res <- lapply(names(assay.lst), function(clstrSet) {
      # for each cluster in the set
      geneCounts <- lapply(assay.lst[[clstrSet]], function(clstr) {
        assay(pb, clstr)
      })
      # sum the multiple clusters in the set
      geneCounts <- Reduce("+", geneCounts)

      # merge with metadata
      #--------------------
      # get names of samples to extract from
      # intersecting between geneExpr and metadata
      ids <- intersect(colnames(geneCounts), rownames(data_constant))
      geneCounts <- geneCounts[, ids, drop = FALSE]
      colnames(geneCounts) <- paste0(clstrSet, "_", colnames(geneCounts))

      # subset and combine metadata
      aggr_means <- as.data.frame(get_metadata_aggr_means(pb))
      include <- aggr_means[[metadata(pb)$agg_pars$by[1]]] %in% assay.lst[[clstrSet]]
      aggr_means <- aggr_means[include, , drop = FALSE]

      lvls <- levels(droplevels(aggr_means[[metadata(pb)$agg_pars$by[1]]]))
      dftmp <- lapply(lvls, function(CT) {
        include <- aggr_means[[metadata(pb)$agg_pars$by[1]]] == CT
        tmp <- aggr_means[include, , drop = FALSE]
        rownames(tmp) <- tmp[[metadata(pb)$agg_pars$by[2]]]
        tmp[, all.vars(formula), drop = FALSE]
      })
      dftmp <- Reduce("+", dftmp) / length(dftmp)

      dftmp[[metadata(pb)$agg_pars$by[1]]] <- clstrSet
      dftmp[[metadata(pb)$agg_pars$by[2]]] <- rownames(dftmp)
      dftmp$cellCluster <- clstrSet
      dftmp$Sample <- rownames(dftmp)

      df <- merge_metadata(data_constant[ids, , drop = FALSE], dftmp, clstrSet, metadata(pb)$agg_pars$by)

      list(geneCounts = geneCounts, df = df)
    })

    countsMatrix <- do.call(cbind, lapply(res, function(x) x$geneCounts))
    data <- do.call(rbind, lapply(res, function(x) x$df))

    ids <- intersect(colnames(countsMatrix), rownames(data))
    rownames(data) <- colnames(countsMatrix)

    # extract cell counts for specified sets
    n.cells <- c(
      rowSums(cellCounts(pb)[, unlist(assay.lst[[1]]), drop = FALSE]),
      rowSums(cellCounts(pb)[, unlist(assay.lst[[2]]), drop = FALSE])
    )
  } else {
    # extract merge pseudobulk counts for specified assays
    res <- lapply(unlist(assay.lst), function(clstr) {
      geneCounts <- assay(pb, clstr)
      colnames(geneCounts) <- paste0(clstr[1], "_", colnames(geneCounts))

      df <- as.data.frame(colData(pb))
      df$cellCluster <- clstr
      df$Sample <- rownames(df)
      rownames(df) <- colnames(geneCounts)

      list(geneCounts = geneCounts, df = df)
    })

    countsMatrix <- do.call(cbind, lapply(res, function(x) x$geneCounts))
    data <- do.call(rbind, lapply(res, function(x) x$df))
    rownames(data) <- colnames(countsMatrix)

    # extract cell counts for specified assays
    n.cells <- c(cellCounts(pb)[, unlist(assay.lst)])
  }

  # create formula to evalute voom and differential expression
  form <- switch(method,
    "random" = {
      update.formula(formula, ~ 0 + cellCluster + (1 | Sample) + .)
    },
    "fixed" = {
      update.formula(formula, ~ 0 + cellCluster + Sample + .)
    },
    "none" = {
      update.formula(formula, ~ 0 + cellCluster + .)
    }
  )

  n.samples <- length(unique(data$Sample))
  n.cellCluster <- length(unique(data$cellCluster))

  # catch error according to errorsAsWarnings
  vobj <- tryCatch(
    {
      # processing counts with voom or log2 CPM
      processOneAssay(countsMatrix, form, data, n.cells,
        min.cells = min.cells,
        min.count = min.count,
        min.samples = min.samples,
        isCounts = isCounts,
        normalize.method = normalize.method,
        quiet = quiet,
        BPPARAM = BPPARAM, ...
      )
    },
    error = function(e) {
      if (errorsAsWarnings) {
        warning(e$message)
      } else {
        stop(e$message)
      }
      NULL
    }
  )

  if (is.null(vobj)) {
    txt <- "No samples passed the filters.\n  Consider looser cutoffs for min.cells, min.count, min.samples"
    if (errorsAsWarnings) {
      warning(txt)
      return(NULL)
    } else {
      stop(txt)
    }
  }

  if (nrow(vobj) < 4) {
    txt <- paste("Only ", nrow(vobj), " samples passed the filters.\n  Consider looser cutoffs for min.cells, min.count, min.samples")

    if (errorsAsWarnings) {
      warning(txt)
      return(NULL)
    } else {
      stop(txt)
    }
  }

  # since vobj contains a subset of cells, also subset the data
  idx <- match(colnames(vobj), rownames(data))
  data <- data[idx, , drop = FALSE]

  # describe samples dropped by filtering
  n.samples2 <- length(unique(data$Sample))
  n.cellCluster2 <- length(unique(data$cellCluster))

  if (!quiet) {
    message("Initial filtering...\n")
    if (n.samples2 - n.samples != 0) {
      message("  Dropped ", (n.samples - n.samples2), "/", n.samples, " samples\n")
    }
    if (n.cellCluster2 - n.cellCluster != 0) {
      message("  Dropped ", (n.cellCluster - n.cellCluster2), "/", n.cellCluster, " cell clusters\n")
    }
  }

  # If paired analysis is requested, and only one example of a Sample is found
  if (method %in% c("random", "fixed")) {
    if (!quiet) {
      message("Filtering for paired samples...\n")
    }

    # drop Samples and cellCluster with only 1 example,
    # 	continue until no more changes are made
    data2 <- data

    n_remaining <- nrow(data2)

    for (i in seq(1, 100)) {
      tab <- table(data2$Sample) > 1
      keep <- data2$Sample %in% names(tab)[tab]
      data2 <- data2[keep, , drop = FALSE]
      tab <- table(data2$cellCluster) > 1
      keep <- data2$cellCluster %in% names(tab)[tab]
      data2 <- data2[keep, , drop = FALSE]

      if (nrow(data2) == n_remaining) break
      n_remaining <- nrow(data2)
    }
    if (i == 100 || n_remaining == 0) {
      txt <- "No samples remain after filtering"

      if (errorsAsWarnings) {
        warning(txt)
        return(NULL)
      } else {
        stop(txt)
      }
    }
    # retain only these samples
    idx <- match(rownames(data2), rownames(data))
    data <- droplevels(data[idx, , drop = FALSE])
    vobj <- vobj[, idx, drop = FALSE]

    if (collapse) {
      # keep only retained cellCluster
      assay.lst$test <- assay.lst$test[assay.lst$test %in% data$cellCluster]
      assay.lst$baseline <- assay.lst$baseline[assay.lst$baseline %in% data$cellCluster]
    }

    # describe samples dropped by filtering	on pairs
    n.samples3 <- length(unique(data$Sample))
    n.cellCluster3 <- length(unique(data$cellCluster))

    if (!quiet) {
      if (n.samples3 - n.samples2 != 0) {
        message("  Dropped ", (n.samples2 - n.samples3), "/", n.samples2, " samples\n")
      }
      if (n.cellCluster3 - n.cellCluster2 != 0) {
        message("  Dropped ", (n.cellCluster2 - n.cellCluster3), "/", n.cellCluster2, " cell clusters\n")
      }
    }

    if (!collapse & min(vapply(assay.lst, length, numeric(1))) == 0) {
      txt <- "Insufficient cellClusters retained after filtering"

      if (errorsAsWarnings) {
        warning(txt)
        return(NULL)
      } else {
        stop(txt)
      }
    }
  }

  if (collapse) {
    L <- makeContrastsDream(form, data, contrasts = contrasts)
  } else {
    # specify contrasts
    test <- paste0("(", paste0("`cellCluster", assay.lst$test, "`", collapse = " + "), ") / ", length(assay.lst$test))
    baseline <- paste0("(", paste0("`cellCluster", assay.lst$baseline, "`", collapse = " + "), ") / ", length(assay.lst$baseline))

    L <- makeContrastsDream(form, data, contrasts = c(compare = paste(test, "-", baseline)))
  }

  # catch error according to errorsAsWarnings
  fit <- tryCatch(
    {
      # perform differential expression regression analysis
      dream(vobj, form, data, L = L, BPPARAM = BPPARAM, ..., quiet = TRUE)
    },
    error = function(e) {
      if (errorsAsWarnings) {
        warning(e$message)
      } else {
        stop(e$message)
      }
      NULL
    }
  )
  if (is.null(fit)) {
    return(NULL)
  }

  # borrow information across genes with the Empirical Bayes step
  fit <- eBayes(fit, robust = robust, trend = !vobj$isCounts)

  # return fit
  fit
}
