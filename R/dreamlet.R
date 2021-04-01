# Gabriel Hoffman
# April 1, 2021
#
# deamlet uses linear mixed models in dream to perform differential expression in single cell data
  
# depends on limma, edgeR, variancePartition

#
# need to adapt voomByDreamWeights to allow per sample weights
# see use of sample weights here: 
# https://rdrr.io/bioc/limma/src/R/voomWithQualityWeights.R



#' Differential expression for each assay
#'
#' Perform differential expression for each assay using linear mixed models
#'
#' @param sceObj SingleCellExperiment object 
#' @param formula regression formula for differential expression analysis
#' @param L contrast matrix specifying a linear combination of fixed effects to test
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param robust logical, use eBayes method that is robust to outlier genes
#' @param method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM: parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import SingleCellExperiment variancePartition BiocParallel limma
#'
#' @export
#'
dreamlet = function( sceObj, formula, L, min.cells = 10, robust=TRUE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	stopifnot( is(sceObj, 'SingleCellExperiment'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	info = as.data.frame(colData(sceObj))

	# for each cell type
	fitList = lapply( assays(sceObj), function(k){
		
		# get data from assay k
		y = assay(sceObj, k)

		# nCells = extract from, sceObj

		# samples to include of they have enough observed cells
		include = nCells >= min.cells
		y = y[,include] 

		# per sample weights based on cell counts in sceObj
		w = nCells[include] #weights by cell types

		# convert vector of sample weights to full matrix
		# each gene is weighted the same
		weights = asMatrixWeights(w, dims=c(nrow(y), ncol(y)))

		if( isCounts ){

			# Get count data and normalize
        	y = suppressMessages(DGEList(y, remove.zeros = TRUE))
        	y = calcNormFactors(y, method = normalize.method )

			# get samples with enough cells
			# filter genes
			design = model.matrix( subbars(formula), info)
			y = filterByExpr(y, design)

			# create EList object storing gene expression and sample weights
			y = new("EList", list(E=y, weights = weights))

			# since the sample weights are already in y, don't need to 
			# explicitly consider them here.
			geneExpr = voomWithDreamWeights( y, formula, info, BPPARAM=BPPARAM,..., save.plot=TRUE)

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

		# fit linear mixed model for each gene
		fit = dream(geneExpr, formula, info, L=L, BPPARAM=BPPARAM,...)

		# borrow information across genes with the Empircal Bayes step
		eBayes(fit, robust=robust, trend=trend)
	})

	# name each result by the assay name
	names(fitList) = assays(sceObj)

	fitList
}










