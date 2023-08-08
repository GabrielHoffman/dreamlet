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
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.signal minimum signal value for a gene to be considered expressed in a sample.  Proper value for this cutoff depends on the type of signal value  
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts for a gene to be
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam}}
#'   object specifying how aggregation should be parallelized.
#' @param verbose logical. Should information on progress be reported?
#'
#' @return a \code{dreamletProcessedData} object 
#' 
#' @details 
#' The \code{dreamlet} workflow can also be applied to non-count data. In this case, a signal is averaged across all cells from a given sample and cell type. Here \code{aggregateNonCountSignal()} performs the roles of \code{aggregateToPseudoBulk()} followed by \code{processAssays()} but using non-count data.
#'
#' For each cell cluster, samples with at least \code{min.cells} are retained. Only clusters with at least \code{min.samples} retained samples are kept. Features are retained if they have at least \code{min.signal} in at least min.prop fraction of the samples.
#'  
#' The precision of a measurement is the inverse of its sampling variance. The precision weights are computed as \code{1/sem^2}, where \code{sem = sd(signal) / sqrt(n)}, \code{signal} stores the values averaged across cells, and \code{n} is the number of cells. 
#' 
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#' 
#' data(example_sce)
#' 
#' # create pseudobulk for each sample and cell cluster
#' # using non-count signal
#' pb.signal <- aggregateNonCountSignal(example_sce, 
#'    assay = "logcounts",    
#'    cluster_id = 'cluster_id', 
#'    sample_id = 'sample_id',
#'    verbose=FALSE)
#' 
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data 
#' res.dl = dreamlet( pb.signal, ~ group_id)
#' @export
#' @importClassesFrom limma EList
#' @importFrom MatrixGenerics rowVars
aggregateNonCountSignal = function(sce, assay = NULL, sample_id = NULL, cluster_id = NULL, min.cells = 10, min.signal = 0.01,  min.samples = 4, min.prop = 0.4, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)){

	# average signal across cells in a sample_id
	pb <- aggregateToPseudoBulk(sce, 
	    assay = assay,    
	    cluster_id = cluster_id, 
	    sample_id = sample_id,
	    verbose=FALSE, 
	    fun = "mean",
	    BPPARAM = BPPARAM,
	    checkValues = FALSE)

	# get the standard error of the mean for each pseudobulk unit
	pb.sem <- aggregateToPseudoBulk(sce, 
	    assay = assay,    
	    cluster_id = cluster_id, 
	    sample_id = sample_id,
	    verbose=FALSE, 
	    fun = "sem",
	    BPPARAM = BPPARAM,
	    checkValues = FALSE)

	# get number of units used for pseudobulk
	# this allows is to distinguish between 
	# a signal of zero due to 
	# 1) no observed cells
	# 2) the mean of multiple cells is zero
	pb.number <- aggregateToPseudoBulk(sce, 
	    assay = assay,    
	    cluster_id = cluster_id, 
	    sample_id = sample_id,
	    verbose=FALSE, 
	    fun = "number",
	    BPPARAM = BPPARAM,
	    checkValues = FALSE)
	
	# extract metadata shared across assays
	data_constant = droplevels(as.data.frame(colData(pb)))

	# Extract signal for each cell type as lists
	resList = lapply(assayNames(pb), function(CT){

		# mean signal for each
		signal = assay(pb, CT)

		# precision weights are in inverse variance. so 1/sem^2
		w = 1 / assay(pb.sem, CT)^2

		# number of cells used for pseudobulk aggregation
		number = assay(pb.number, CT)

		# create EList with signal and precision weights
		obj = new('EList', list(E = as.matrix(signal), weights = as.matrix(w)))

		# inclusion criteria for samples
		include = (number[1,] >= min.cells)
		obj = obj[,include,drop=FALSE]

		# inclusion criteria for genes
		keep = apply(obj$E, 1, function(x){
			sum(x >= min.signal) >= min.prop*length(x)
			} )
		obj = obj[keep,,drop=FALSE]

		# if weight is Inf, set to largest finite value per row
		if( any(!is.finite(obj$weights)) ){
			for(i in seq(nrow(obj$weights))){
				if( any(!is.finite(obj$weights[i,])) ){
					idx = which(!is.finite(obj$weights[i,]))
					if( length(idx) > 0){
						obj$weights[i,idx] = max(obj$weights[i, -idx])
					}
				}
			}
		}

		# if there are too few remaining samples
		if( ncol(obj) < min.samples ){
			return( NULL )
		}

		obj
	})
	names(resList) = assayNames(pb)

	# remove empty assays
	resList = resList[!vapply(resList, is.null, FUN.VALUE=logical(1))]

	# return signal as dreamletProcessedData object
	new("dreamletProcessedData", 
		resList, 
		data = data_constant, 
		metadata = metadata(pb)$aggr_means, 
		by = metadata(pb)$agg_pars$by)
}

