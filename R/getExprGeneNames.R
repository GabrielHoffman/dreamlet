

#' Get list of expressed genes for each assay
#'
#' Get list of expressed genes for each assay using same filters as \code{processAssays()}.
#'
#' @param sceObj SingleCellExperiment object
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(sceObj)}
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count minimum number of reads for a gene to be considered expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts for a gene to be retained
#' @param min.total.count minimum total count required per gene for inclusion
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
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
#' @export
getExprGeneNames <- function(sceObj, assays = assayNames(sceObj), min.cells = 5, min.count = 5, min.samples = 4, min.prop = .4, min.total.count = 15, normalize.method = "TMM"){

  # checks
  stopifnot(is(sceObj, "SingleCellExperiment"))

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
    startTime <- Sys.time()

    y <- as.matrix(assay(sceObj, k))

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

    y = y[, rownames(data), drop = FALSE]

  	# samples to include of they have enough observed cells
  	include <- (n.cells >= min.cells)

  	# if no samples are retained
  	if (sum(include) == 0) {
  		return(NULL)
  	}

  	# subset expression and data
  	y <- y[, include, drop = FALSE]
  	data <- droplevels(data[include, , drop = FALSE])

  	# if there are too few remaining samples
  	if (nrow(data) < min.samples | nrow(y) == 0) {
  		return(NULL)
  	}

  	# Get count data and normalize
  	y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
  	y <- calcNormFactors(y, method = normalize.method)

  	# get samples with enough cells
  	# filter genes
  	# design: model.matrix( subbars(formula), data)
  	# Design often includes batch and donor, which are very small
  	#   this causes too many genes to be retained
  	keep <- suppressWarnings(
    filterByExpr(y, 
      min.count = min.count, 
      min.prop = min.prop, 
      min.total.count = min.total.count)
    )
  	names(keep)[keep]
	})
	names(resList) <- assays
	resList
}



