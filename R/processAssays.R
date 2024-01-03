#' Processing expression data from assay
#'
#' For raw counts, filter genes and samples, then estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param y matrix of counts or log2 CPM
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param n.cells array of cell count for each sample
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count minimum number of reads for a gene to be considered expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param span Lowess smoothing parameter using by \code{variancePartition::voomWithDreamWeights()}
#' @param quiet show messages
#' @param weights matrix of precision weights 
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @return \code{EList} object storing log2 CPM and precision weights
#'
#' @seealso \code{processAssays()}
#' @importFrom BiocParallel SerialParam
#' @importClassesFrom limma EList
#' @importFrom variancePartition voomWithDreamWeights
#' @importFrom edgeR calcNormFactors filterByExpr DGEList
#' @importFrom methods is new
#' @importFrom stats model.matrix var
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @importFrom lme4 subbars
#' @importFrom MatrixGenerics colMeans2
#'
processOneAssay <- function(y, formula, data, n.cells, min.cells = 5, min.count = 5, min.samples = 4, min.prop = .4, isCounts = TRUE, normalize.method = "TMM", span = "auto", quiet = TRUE, weights=NULL, BPPARAM = SerialParam(), ...) {

  checkFormula(formula, data)

  if (is.null(n.cells)) {
    stop("n_cells must not be NULL")
  }
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  # samples to include of they have enough observed cells
  include <- (n.cells >= min.cells)

  # if no samples are retained
  if (sum(include) == 0) {
    return(NULL)
  }

  # subset expression and data
  y <- y[, include, drop = FALSE]
  data <- droplevels(data[include, , drop = FALSE])
  # n.cells <- n.cells[include]

  # if there are too few remaining samples
  if (nrow(data) < min.samples | nrow(y) == 0) {
    return(NULL)
  }

  if ( ! isCounts ) {
    stop("isCounts = FALSE is not currently supported")
  }

  # Get count data and normalize
  y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
  y <- calcNormFactors(y, method = normalize.method)

  # drop any constant terms from the formula
  formula <- removeConstantTerms(formula, data)

  # Drop variables in a redundant pair
  formula <- dropRedundantTerms(formula, data)

  # get samples with enough cells
  # filter genes
  # design: model.matrix( subbars(formula), data)
  # Design often includes batch and donor, which are very small
  #   this causes too many genes to be retained
  keep <- suppressWarnings(filterByExpr(y, min.count = min.count, min.prop = min.prop))

  # sample-level weights based on cell counts and mean library size
  if ( ! is.null(weights) & ! is(weights, "function") ) {
    precWeights <- weights[rownames(y)[keep],colnames(y)]
  } else {
    precWeights <- rep(1, ncol(y))
  }

  # if no genes are kept
  if( sum(keep) == 0) return(NULL)

  geneExpr <- voomWithDreamWeights(y[keep, ], formula, data, weights = precWeights, BPPARAM = BPPARAM, ..., save.plot = TRUE, quiet = quiet, span = span, hideErrorsInBackend = TRUE)

  # if no genes are succeed
  if( nrow(geneExpr) == 0) return(NULL)

  # save formula used after dropping constant terms
  if (!is.null(geneExpr)) geneExpr$formula <- formula
  if (!is.null(geneExpr)) geneExpr$isCounts <- isCounts

  geneExpr
}







#' Processing SingleCellExperiment to dreamletProcessedData
#'
#' For raw counts, estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells.
#'
#' @param sceObj SingleCellExperiment object
#' @param formula regression formula for differential expression analysis
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(sceObj)}
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count minimum number of reads for a gene to be considered expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts for a gene to be retained
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param span Lowess smoothing parameter using by \code{variancePartition::voomWithDreamWeights()}
#' @param quiet show messages
#' @param weightsList list storing matrix of precision weights for each cell type. If \code{NULL} precision weights are set to 1
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @return Object of class \code{dreamletProcessedData} storing voom-style normalized expression data
#'
#' @details  For each cell cluster, samples with at least \code{min.cells} are retained. Only clusters with at least \code{min.samples} retained samples are kept. Genes are retained if they have at least \code{min.count} reads in at least \code{min.prop} fraction of the samples.  Current values are reasonable defaults, since genes that don't pass these cutoffs are very underpowered for differential expression analysis and only increase the multiple testing burden.  But values of \code{min.cells = 2} and \code{min.count = 2} are also reasonable to include more genes in the analysis.
#'
#' The precision weights are estimated using the residuals fit from the specified formula.  These weights are robust to changes in the formula as long as the major variables explaining the highest fraction of the variance are included.  
#' 
#' If \code{weightsList} is \code{NULL}, precision weights are set to 1 internally.
#'
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
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~ group_id)
#'
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data
#' res.dl <- dreamlet(res.proc, ~ group_id)
#
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors metadata as.data.frame
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays assay
#
#' @export
processAssays <- function(sceObj, formula, assays = assayNames(sceObj), min.cells = 5, min.count = 5, min.samples = 4, min.prop = .4, isCounts = TRUE, normalize.method = "TMM", span = "auto", quiet = FALSE, weightsList = NULL, BPPARAM = SerialParam(), ...) {
  # checks
  stopifnot(is(sceObj, "SingleCellExperiment"))
  stopifnot(is(formula, "formula"))

  # check colnames of SCE
  if (is.null(colnames(sceObj))) {
    stop("colnames(sceObj) is NULL.  Column names are needed for internal filtering")
  }

  # extract metadata shared across assays
  data_constant <- droplevels(as.data.frame(colData(sceObj)))

  # check if assays are valid
  if (any(!assays %in% assayNames(sceObj))) {
    idx <- which(!assays %in% assayNames(sceObj))
    txt <- paste("Assays are not found in dataset:", paste(head(assays[idx]), collapse = ", "))
    stop(txt)
  }

  # extract cell counts
  n.cells_full <- cellCounts(sceObj)

  # extract all unique colnames
  colNamesAll <- unique(unlist(lapply(assayNames(sceObj), function(x) colnames(assay(sceObj, x)))))

  # check for colnames missing cell counts
  if (any(!colNamesAll %in% rownames(n.cells_full))) {
    stop("Cell counts could not be extracted.\n  Check that colnames(sceObj) or rownames(colData(sceObj))\n  have not been modified manually after running aggregateToPseudoBulk()")
  }

  # for each assay
  resList <- lapply(assays, function(k) {
    if (!quiet) message("  ", k, "...", appendLF = FALSE)
    startTime <- Sys.time()

    y <- assay(sceObj, k)

    # cell counts
    n.cells <- n.cells_full[colnames(y), k, drop = FALSE]

    # merge data_constant (data constant for all cell types)
    # with metadata(sceObj)$aggr_means (data that varies)
    data <- merge_metadata(
      data_constant,
      get_metadata_aggr_means(sceObj),
      k,
      metadata(sceObj)$agg_pars$by
    )

    if( ! is.null(weightsList) ){
      # use specified weights
      weights <- weightsList[[k]][,rownames(data),drop=FALSE]
    }else{
      # use cell count weights
      # weights <- n.cells[rownames(data),]
      # weights are 1
      weights <- matrix(1, nrow = nrow(y), ncol = nrow(data), byrow=TRUE)
      colnames(weights) <- rownames(data)
      rownames(weights) <- rownames(y)
    }

    # processing counts with voom or log2 CPM
    res <- processOneAssay(y[, rownames(data), drop = FALSE],
      formula = formula,
      data = data,
      n.cells = n.cells[rownames(data), , drop = FALSE],
      min.cells = min.cells,
      min.count = min.count,
      min.samples = min.samples,
      min.prop = min.prop,
      isCounts = isCounts,
      normalize.method = normalize.method,
      span = span,
      weights = weights,
      BPPARAM = BPPARAM, 
      ...
    )

    if (!quiet) message(format(Sys.time() - startTime, digits = 2))

    res
  })
  names(resList) <- assays

  exclude <- vapply(resList, is.null, FUN.VALUE = logical(1))

  if (length(names(resList)[exclude]) > 0) {
    warning("Not enough samples retained or model fit fails: ", paste(names(resList)[exclude], collapse = ", "))
  }

  # remove empty assays
  resList <- resList[!exclude]

  # Handle errors
  #--------------

  # get error messages
  error.initial <- lapply(resList, function(x) {
    x$error.initial
  })
  names(error.initial) <- names(resList)
  errors <- lapply(resList, function(x) {
    attr(x, "errors")
  })
  names(errors) <- names(resList)

  # extract details
  df_details <- lapply(names(resList), function(id) {
    data.frame(
      assay = id,
      n_retain = ncol(resList[[id]]),
      formula = Reduce(paste, deparse(resList[[id]]$formula)),
      formDropsTerms = !equalFormulas(resList[[id]]$formula, formula),
      n_genes = nrow(resList[[id]]),
      n_errors = length(attr(resList[[id]], "errors")),
      error_initial = ifelse(is.null(resList[[id]]$error.initial), FALSE, TRUE)
    )
  })
  df_details <- do.call(rbind, df_details)

  ndrop <- sum(df_details$formDropsTerms)

  if (ndrop > 0) {
    warning("Terms dropped from formulas for ", ndrop, " assays.\n Run details() on result for more information")
  }

  failure_frac <- sum(df_details$n_errors) / sum(df_details$n_genes)

  if( is.nan(failure_frac) ){
    stop("All models failed.  Consider changing formula")
  }

  if( failure_frac > 0 ){
    txt <- paste0("\nOf ", format(sum(df_details$n_genes), big.mark=','), " models fit across all assays, ", format(failure_frac*100, digits=3), "% failed\n")
    message(txt)
  }

  new("dreamletProcessedData",
    resList,
    data = data_constant,
    # metadata = metadata(sceObj)$aggr_means,
    metadata = get_metadata_aggr_means(sceObj),    
    by = metadata(sceObj)$agg_pars$by,
    df_details = df_details,
    errors = errors,
    error.initial = error.initial
  )
}

# merge data_constant (data constant for all cell types)
# with metadata(sceObj)$aggr_means (data that varies)
#' @importFrom dplyr filter_at
merge_metadata <- function(dataIn, md, cellType, by) {
  # PASS R CMD check
  cell <- NULL

  data <- merge(dataIn,
    dplyr::filter_at(md, by[1], ~ . == cellType),
    by.x = "row.names",
    by.y = by[2]
  )
  rownames(data) <- data$Row.names
  id <- rownames(dataIn)[rownames(dataIn) %in% rownames(data)]
  data <- data[id, , drop = FALSE]
  droplevels(data)
}