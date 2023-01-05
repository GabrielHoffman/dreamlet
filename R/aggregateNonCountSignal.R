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
#' @details For standard analysis of count data, \code{aggregateToPseudoBulk()} create pseudobulk counts and \code{processAssays()} performs folder and estimates precision weights.  Yet the \code{dreamlet} workflow can also be applied to non-count data. In this case, a signal is averaged across all cells from a given sample and cell type.
#'  
#' The precision of a measurement is the inverse of its sampling variance. The precision weights are computed as \code{1/sem^2}, where \code{sem = sd(signal) / sqrt(n)}, \code{signal} stores the values averaged across cells, and \code{n} is the number of cells.  The weights are modified from this value in 3 edge cases.
#'
#' 1) if \code{n==0} so that no cells are observed for a given cell type and sample, then the signal is set to NA, and weight is set to 0.
#'
#' 2) if \code{n==1}, then there is only one cell and the sem is zero. The weight is set to match the lowest non-zero value for that feature
#'
#' 3) if \code{weight == Inf}, then sem = 0 so set weight is set to match the largest non-zero value for that feature
#' 
#' @export
#' @import limma 
#' @importFrom MatrixGenerics rowVars
aggregateNonCountSignal = function(sce, assay = NULL, sample_id = NULL, cluster_id = NULL, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)){

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
	# 2) the mean of multiple cells iz erpo
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
	pmetadata = data.frame()
	pkeys = array()

	# Extract signal as lists
	resList = lapply(assayNames(pb), function(CT){

		# mean signal for each
		signal = assay(pb, CT)

		# precision weights are in inverse variance. so 1/sem^2
		w = 1 / assay(pb.sem, CT)^2

		# number of cells used for pseudobulk aggregation
		number = assay(pb.number, CT)

		# if zero counts, set E = NA, w = NA
		idx = which(number == 0)
		if( length(idx) > 0){
			signal[idx] = NA
			w[idx] = 0
		}

		# if 1 count, sem is currently set to NA.  
		# Want E = signal, w = lowest weight observed value per row
		idx = which(number[1,] == 1)
		if( length(idx) > 0){
			w[,idx] = apply(w, 1, function(b){
				b = b[b>0 & is.finite(b)]
				ifelse(length(b) > 0, min(b), 1)
				})
		}

		# if weight is Inf, set to largest finite value per row
		if( any(!is.finite(w)) ){
			for(i in seq(nrow(w))){
				if( all(!is.finite(w[i,])) ){
					w[i,] = 0
				}else{ 
					idx = which(!is.finite(w[i,]))
					if( length(idx) > 0){
						w[i,idx] = max(w[i, -idx])
					}
				}
			}
		}
		
		# create EList with signal and precision weights
		obj = new('EList', list(E = as.matrix(signal), weights = as.matrix(w)))

		# keep rows that satisfy
		# variance of signal > 0, and maximum weight > 0
		keep = (rowVars(obj$E, na.rm=TRUE) > 0) & (apply(obj$w, 1, max) > 0) 

		obj[keep,]
		})
	names(resList) = assayNames(pb)

	# return signal as dreamletProcessedData object
	new("dreamletProcessedData", resList, data = data_constant, metadata = pmetadata, pkeys=pkeys)
}






# ids = c('Rush-35', 'Rush-38', 'Rush-39')

# assay(pb, CT)['ADNP2(+)', ids]
# assay(pb.sem, CT)['ADNP2(+)',ids]
# assay(pb.number, CT)['ADNP2(+)',ids]




# w[1:4, ids]



