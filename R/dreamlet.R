# Gabriel Hoffman
# April 1, 2021
#
# deamlet uses linear mixed models in dream to perform differential expression in single cell data
  
# depends on limma, edgeR, variancePartition

#
# need to adapt voomByDreamWeights to allow per sample weights
# see use of sample weights here: 
# https://rdrr.io/bioc/limma/src/R/voomWithQualityWeights.R


# 1) 
# resList = dreamlet( sceObj)

# 2)
# # process data
# # dreamletProcessedData
# procData = processAssays( sceObj)

# # dreamlet on each processed assay
# fitList = dreamlet( procData)

# fitVarPart( fitList$data)
# fitVarPart( procData)

# zenith( fitList$fit )


# processAssay() needs to take nCells as arguments

# should dreamletProcessedData store a SingleCellExperiment?



#' Class dreamletProcessedData
#'
#' Class \code{dreamletProcessedData} 
#'
#' @name dreamletProcessedData-class
#' @rdname dreamletProcessedData-class
#' @exportClass dreamletProcessedData
setClass("dreamletProcessedData", contains="list", slots = c(data = 'data.frame'))
	# representation("list", data="data.frame"))







#' Processing expression data from assay
#'
#' For raw counts, estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param y matrix of counts or log2 CPM
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param n.cells array of cell count for each sample
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel 
#' @import limma 
#' @importFrom variancePartition voomWithDreamWeights
#' @importFrom edgeR calcNormFactors filterByExpr DGEList 
#' @importFrom lme4 subbars  
#' @importFrom methods is new
#' @importFrom stats model.matrix var
#' @importFrom SummarizedExperiment as.data.frame colData assays
#'
#' @export
processOneAssay = function( y, formula, data, n.cells, min.cells = 10, isCounts = TRUE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# nCells = extract from y

	# samples to include of they have enough observed cells
	include = (n.cells >= min.cells)

	# if no samples are retained
	if( sum(include) == 0){
		return(NULL)
	}

	y = y[,include] 

	# per sample weights based on cell counts in sceObj
	w = n.cells[include] #weights by cell types

	# convert vector of sample weights to full matrix
	# each gene is weighted the same
	weights = asMatrixWeights(w, dim=c(nrow(y), ncol(y)))

	if( isCounts ){

		# Get count data and normalize
    	y = suppressMessages(DGEList(y, remove.zeros = TRUE))
    	y = calcNormFactors(y, method = normalize.method )

		# get samples with enough cells
		# filter genes
		design = model.matrix( subbars(formula), data)
		y = filterByExpr(y, design)

		# create EList object storing gene expression and sample weights
		y = new("EList", list(E=y, weights = weights))

		# since the sample weights are already in y, don't need to 
		# explicitly consider them here.
		geneExpr = voomWithDreamWeights( y, formula, data, BPPARAM=BPPARAM,..., save.plot=TRUE)

		# voom dots and curves are saved here
		# geneExpr$voom.line
		# geneExpr$voom.xy

		trend = FALSE
	}else{

		# only include genes that show variation
		include = apply(y, 1, var) > 0

		# if data is already log2 CPM
		# create EList object storing gene expression and sample weights
		geneExpr = new("EList", list(E=y[include,], weights = weights[include,]))

		# since precision weights are not used, use the trend in the eBayes step
		trend = TRUE
	}

	list(geneExpr = geneExpr, trend=trend)
}



#' Processing SingleCellExperiment to dreamletProcessedData
#'
#' For raw counts, estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param sceObj SingleCellExperiment object 
#' @param formula regression formula for differential expression analysis
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel  
#' @importFrom SummarizedExperiment as.data.frame colData assays assay
#' @importFrom S4Vectors metadata
#'
#' @export
processAssays = function( sceObj, formula, min.cells = 10, isCounts=TRUE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# checks
	stopifnot( is(sceObj, 'SingleCellExperiment'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data = as.data.frame(colData(sceObj))

	# for each assay
	resList = lapply( names(assays(sceObj)), function(k){

		y = assay(sceObj, k)
		n.cells = metadata(sceObj)$n_cells[k,colnames(y)]

		# processing counts with voom or log2 CPM
		processOneAssay(y, formula, data, n.cells, min.cells, isCounts, normalize.method, BPPARAM=BPPARAM,...)
	})
	names(resList) = names(assays(sceObj))

	new("dreamletProcessedData", resList, data=data)
}



#' Differential expression for each assay
#'
#' Perform differential expression for each assay using linear mixed models
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param L contrast matrix specifying a linear combination of fixed effects to test
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param isCounts logical, indicating if data is raw counts
#' @param robust logical, use eBayes method that is robust to outlier genes
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel  
#' @importFrom SummarizedExperiment as.data.frame colData assays
#'
#' @export
setGeneric("dreamlet", 
	function( x, formula, data, L, min.cells = 10, isCounts=TRUE, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	standardGeneric("dreamlet")
})




#' @importFrom SummarizedExperiment as.data.frame colData assays assay
#' @importFrom variancePartition dream eBayes
#' @export
#' @rdname dreamlet
#' @aliases dreamlet,SingleCellExperiment-method
setMethod("dreamlet", "SingleCellExperiment",
	function( x, formula, data, L, min.cells = 10, isCounts=TRUE, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'SingleCellExperiment'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data = as.data.frame(colData(x))

	# for each assay
	resList = lapply( names(assays(x)), function(k){

		message("\rAssay: ", k)
		
		# get data from assay k
		y = assay(x, k)
		n.cells = metadata(x)$n_cells[k,colnames(y)]

		# processing counts with voom or log2 CPM
		res = processOneAssay(y, formula, data, n.cells, min.cells, isCounts, normalize.method, BPPARAM=BPPARAM,...)

		# if samples are retained after filtering
		if( ! is.null(res) ){

			# fit linear (mixed) model for each gene
			# only include samples from data that are retained in res$geneExpr
			# TODO include L now
			fit = dream( res$geneExpr, formula, data[colnames(res$geneExpr),], BPPARAM=BPPARAM, quiet=TRUE,...)

			# if model is degenerate
			if( ! any(is.na(fit$sigma)) ){
				# borrow information across genes with the Empircal Bayes step
				fit = eBayes(fit, robust=robust, trend=res$trend)
			}else{	
				fit = NULL
			}
		}else{
			fit = NULL
		}

		list(fit = fit, data = res)
	})
	# name each result by the assay name
	names(resList) = names(assays(x))

	# create list of all fit objects
	fitList = lapply(resList, function(obj) obj$fit)
	names(fitList) = names(assays(x))

	# create list of all data objects
	dataList = lapply(resList, function(obj) obj$data)
	names(dataList) = names(assays(x))

	list(fit = fitList, data = new("dreamletProcessedData", dataList, data=data) )
})






#' @importFrom SummarizedExperiment as.data.frame colData assays
#' @export
#' @rdname dreamlet
#' @aliases dreamlet,dreamletProcessedData-method
setMethod("dreamlet", "dreamletProcessedData",
	function( x, formula, data, L, min.cells = 10, isCounts=TRUE, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	if( missing(data) ){
		data = as.data.frame(colData(x))
	}

	# for each assay
	resList = lapply( x, function( procData ){

		# get names of samples to extract from metadata
		ids = colnames(procData$geneExpr)

		# fit linear mixed model for each gene
		# TODO add , L=L
		fit = dream( procData$geneExpr, formula, data[ids,], BPPARAM=BPPARAM,...)

		# if model is degenerate
		if( ! any(is.na(fit$sigma)) ){
			# borrow information across genes with the Empircal Bayes step
			fit = eBayes(fit, robust=robust, trend=procData$trend)
		}else{	
			fit = NULL
		}
		fit
	})
	# name each result by the assay name
	names(resList) = names(x)

	resList
})




# setGeneric("colData",
# 	function(x,...){		
# 	standardGeneric("colData")
# })


#' Extract colData from dreamletProcessedData
#' 
#' Extract colData from dreamletProcessedData
#' @param x A dreamletProcessedData object
#' @param ... other arguments
# @export
setMethod("colData", "dreamletProcessedData",
	function(x,...){
		x@data
})



#' Variance Partition analysis for each assay
#'
#' Perform Variance Partition analysis  for each assay
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel  
#'
#' @export
setGeneric("fitVarPart", 
	function( x, formula, data, BPPARAM = bpparam(),...){

	standardGeneric("fitVarPart")
})




#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom SummarizedExperiment as.data.frame colData assays
#' @export
#' @rdname fitVarPart
#' @aliases fitVarPart,dreamletProcessedData-method
setMethod("fitVarPart", "dreamletProcessedData",
	function( x, formula, data, BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))
	
	# extract metadata shared across assays
	if( missing(data) ){
		data = as.data.frame(colData(x))
	}

	# for each assay
	resList = lapply( x, function( procData ){

		# get names of samples to extract from metadata
		ids = colnames(procData$geneExpr)

		# fit linear mixed model for each gene
		# TODO add , L=L
		fitExtractVarPartModel( procData$geneExpr, formula, data[ids,], BPPARAM=BPPARAM,...)
	})
	# name each result by the assay name
	names(resList) = names(x)

	resList
})














