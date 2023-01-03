# Gabriel Hoffman
# Jan 3, 2023



#' Aggregation of single-cell signals
#' 
#' Aggregation of single-cell to pseudobulk data for non-count data.  
#' 
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying the assay slot to use as 
#'   input data. Defaults to the 1st available (\code{assayNames(x)[1]}).
# @param by character vector specifying which 
#   \code{colData(x)} columns to summarize by (at most 2!).
#' @param sample_id character string specifying which variable to use as sample id
#' @param cluster_id character string specifying which variable to use as cluster id
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam}}
#'   object specifying how aggregation should be parallelized.
#' @param verbose logical. Should information on progress be reported?
#'
#' @return a \code{dreamletProcessedData} object 
#' 
#' @details For standard analysis of count data, \code{aggregateToPseudoBulk()} create pseudobulk counts and \code{processAssays()} performs folder and estimates precision weights.  Yet the \code{dreamlet} workflow can also be applied to non-count data. In this case, a signal is averaged across all cells from a given sample and cell type.  And the precision weights step is set skipped.          
#' 
#' @export
aggregateNonCountSignal = function(sce, assay = NULL, sample_id = NULL, cluster_id = NULL, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)){

	# average signal across cells in a sample_id
	pb <- aggregateToPseudoBulk(sce, 
	    assay = assay,    
	    cluster_id = cluster_id, 
	    sample_id = sample_id,
	    verbose=FALSE, 
	    fun = "mean",
	    BPPARAM = BPPARAM)

	# extract metadata shared across assays
	data_constant = droplevels(as.data.frame(colData(pb)))
	pmetadata = data.frame()
	pkeys = array()

	# Extract signal as lists
	resList = lapply(assayNames(pb), function(CT){
		assay(pb, CT)
		})
	names(resList) = assayNames(pb)

	# return signal as dreamletProcessedData object
	new("dreamletProcessedData", resList, data = data_constant, metadata = pmetadata, pkeys=pkeys)
}
