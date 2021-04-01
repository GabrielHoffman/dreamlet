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




#' Class dreamletProcessedData
#'
#' Class \code{dreamletProcessedData} 
#'
#' @name dreamletProcessedData-class
#' @rdname dreamletProcessedData-class
#' @exportClass dreamletProcessedData
setClass("dreamletProcessedData", representation("list"))




#' Processing expression data from assay
#'
#' For raw counts, estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param y matrix of counts or log2 CPM
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel 
#' @import limma 
#' @import variancePartition
#' @importFrom edgeR calcNormFactors filterByExpr DGEList 
#' @importFrom lme4 subbars  
#' @importFrom methods is new
#' @importFrom stats model.matrix
#'
#' @export
processOneAssay = function( y, formula, data, min.cells = 10, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# nCells = extract from y

	# samples to include of they have enough observed cells
	include = nCells >= min.cells
	y = y[,include] 

	# per sample weights based on cell counts in sceObj
	w = nCells[include] #weights by cell types

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

		# if data is already log2 CPM
		# create EList object storing gene expression and sample weights
		geneExpr = new("EList", list(E=y, weights = weights))

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
#' @param data metadata used in regression formula
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel SingleCellExperiment  SummarizedExperiment
#'
#' @export
processAssays = function( sceObj, formula, min.cells = 10, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# checks
	stopifnot( is(sceObj, 'SingleCellExperiment'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data = as.data.frame(colData(sceObj))

	# for each assay
	resList = lapply( assays(sceObj), function(k){

		y = assay(sceObj, k)

		# processing counts with voom or log2 CPM
		processOneAssay(y, formula, data, min.cells, normalize.method, BPPARAM=BPPARAM,...)
	})

	names(resList) = names(sceObj)

	new("dreamletProcessedData", resList)
}



#' Differential expression for each assay
#'
#' Perform differential expression for each assay using linear mixed models
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param L contrast matrix specifying a linear combination of fixed effects to test
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param robust logical, use eBayes method that is robust to outlier genes
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel SingleCellExperiment 
#' @import SummarizedExperiment
#'
#' @export
setGeneric("dreamlet", 
	function( x, formula, data, L, min.cells = 10, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	standardGeneric("dreamlet")
})




#' @export
# @rdname dreamlet
# @aliases dreamlet,SingleCellExperiment-method
setMethod("dreamlet", "SingleCellExperiment",
	function( x, formula, data, L, min.cells = 10, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'SingleCellExperiment'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data = as.data.frame(colData(x))

	# for each assay
	resList = lapply( assays(x), function(k){
		
		# get data from assay k
		y = assay(x, k)

		# processing counts with voom or log2 CPM
		res = processOneAssay(y, formula, data, min.cells, normalize.method, BPPARAM=BPPARAM,...)

		# fit linear mixed model for each gene
		fit = dream( res$geneExpr, formula, data, L=L, BPPARAM=BPPARAM,...)

		# borrow information across genes with the Empircal Bayes step
		fit = eBayes(fit, robust=robust, trend=res$trend)

		list(fit = fit, data = res$geneExpr)
	})

	# name each result by the assay name
	names(resList) = assays(x)

	# create list of all fit objects
	fitList = lapply(resList, function(obj) obj$fit)
	names(fitList) = assays(x)

	# create list of all data objects
	dataList = lapply(resList, function(obj) obj$data)
	names(dataList) = assays(x)

	list(fit = fitList, data = new("dreamletProcessedData", dataList) )
})






#' @export
# @rdname dreamlet
# @aliases dreamlet,dreamletProcessedData-method
setMethod("dreamlet", "dreamletProcessedData",
	function( x, formula, data, L, min.cells = 10, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	if( missing(data) ){
		data = as.data.frame(colData(x))
	}

	# for each assay
	resList = lapply( x, function( procData ){

		# fit linear mixed model for each gene
		fit = dream( procData, formula, data, L=L, BPPARAM=BPPARAM,...)

		# borrow information across genes with the Empircal Bayes step
		fit = eBayes(fit, robust=robust, trend=procData$trend)

		list(fit=fit)
	})

	# name each result by the assay name
	names(resList) = assays(x)

	resList
})



#' Variance Partition analysis for each assay
#'
#' Perform Variance Partition analysis  for each assay
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel  
#' @import SummarizedExperiment
#'
#' @export
setGeneric("fitVarPart", 
	function( x, formula, data, BPPARAM = bpparam(),...){

	standardGeneric("fitVarPart")
})




#' @importFrom variancePartition fitExtractVarPartModel
#' @export
setMethod("fitVarPart", "dreamletProcessedData",
	function( x, formula, data, BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))
	stopifnot( is(data, 'data.frame'))

	# for each assay
	resList = lapply( x, function( procData ){

		fitExtractVarPartModel( procData, form, data, BPPARAM = bpparam(), ... )
	})

	# name each result by the assay name
	names(resList) = assays(x)

	resList
})














