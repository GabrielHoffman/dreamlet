# Gabriel Hoffman
# April 1, 2021
#
# deamlet uses linear mixed models in dream to perform differential expression in single cell data
  
# depends on limma, edgeR, variancePartition

#' Differential expression for each assay
#'
#' Perform differential expression for each assay using linear mixed models
#'
#' @param sceObj SingleCellExperiment object 
#' @param formula regression formula for differential expression analysis
#'
#' @import SingleCellExperiment variancePartition
#'
#' @export
#'
dreamlet = function( sceObj, formula){

	# for each cell type
	fitList = lapply( assays(sceObj), function(cellType){
		
		# get data for this cell type
		geneExpr = assay(snObj, cellType)

		# get samples with enough cells
		# filter genes
		geneExpr = filterByExpr(geneExpr)

		# per sample weights
		w = weights by cell types

		# need to adapt voomByDreamWeights to allow per sample weights
		vobj = voomByDreamWeights( geneExpr, form, info, weights=w)

		# need to adapt dream to allow per sample weights
		fit = dream(vobj, form, info, weights=w)

		# finish eBayes
		# only use trend if precion weights are not available
		fit = eBayes(fit, robust, trend)

		# return full results like muscat does
		fit
	})

	# Add ashr on the results of topTable

	fitList
}

