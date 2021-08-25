
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
	function( x, formula, data = colData(x), BPPARAM = SerialParam(),...){

	standardGeneric("fitVarPart")
})




#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom SummarizedExperiment colData assays
#' @importFrom data.table data.table
#' @importFrom S4Vectors DataFrame as.data.frame
#' @export
#' @rdname fitVarPart
#' @aliases fitVarPart,dreamletProcessedData-method
setMethod("fitVarPart", "dreamletProcessedData",
	function( x, formula, data = colData(x), BPPARAM = SerialParam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))
	
	# extract metadata shared across assays
	data_constant = as.data.frame(data)

	# for each assay
	resList = lapply( assayNames(x), function( k ){

		geneExpr = assay(x, k)

		# get names of samples to extract from metadata
		ids = colnames(geneExpr)

		# merge data_constant and pmetadata based on pkeys and assay k
		data2 = merge_metadata(data_constant[ids,,drop=FALSE], metadata(x), pkeys, k)
		data2 = droplevels(data2)

		# drop any constant terms from the formula
		form_mod = removeConstantTerms(formula, data2)

		# fit linear mixed model for each gene
		# TODO add , L=L
		fitExtractVarPartModel( geneExpr, form_mod, data2, BPPARAM=BPPARAM,...,quiet=TRUE)
	})
	# name each result by the assay name
	names(resList) = names(x)

	# Convert results to DataFrame in vpDF
	vplst = lapply( names(resList), function(id){
		data.table(assay = id, gene = rownames(resList[[id]]), data.frame(resList[[id]]))
	})
	df = do.call(rbind, vplst)
	`:=` = NULL # Pass R CMD check
	df[,assay:=factor(assay, names(resList))]
	new("vpDF", DataFrame(df))
})



