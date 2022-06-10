

#' Variance Partition analysis for each assay
#'
#' Perform Variance Partition analysis  for each assay
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(x)}
#' @param quiet show messages
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @return Object of class \code{vpDF} inheriting from \code{DataFrame} storing the variance fractions for each gene and cell type.
#'
#' @examples
#'  
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce, 
#'    assay = "counts",    
#'    cluster_id = 'cluster_id', 
#'    sample_id = 'sample_id',
#'    verbose=FALSE)
#'
#' # voom-style normalization
#' res.proc = processAssays( pb, ~ group_id)
#' 
#' # variance partitioning analysis
#' vp = fitVarPart( res.proc, ~ group_id)
#' 
#' # Show variance fractions at the gene-level for each cell type
#' genes = vp$gene[2:4]
#' plotPercentBars(vp[vp$gene %in% genes,])
#' 
#' # Summarize variance fractions genome-wide for each cell type
#' plotVarPart(vp)
#'
#' @import BiocParallel  
#' @export
setGeneric("fitVarPart", 
	function( x, formula, data = colData(x), assays = assayNames(x), quiet = FALSE, BPPARAM = SerialParam(),...){

	standardGeneric("fitVarPart")
})


# local definition so methods in this file have this class
setClass("dreamletProcessedData", contains="list", slots = c(data = 'data.frame', metadata='data.frame', pkeys="vector"))

#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors DataFrame as.data.frame
#' @importFrom gtools smartbind
#' @export
#' @rdname fitVarPart
#' @aliases fitVarPart,dreamletProcessedData-method
setMethod("fitVarPart", "dreamletProcessedData",
	function( x, formula, data = colData(x), assays = assayNames(x), quiet = FALSE, BPPARAM = SerialParam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))

	# check if assays are valid
	if( any( ! assays %in% assayNames(x)) ){
		idx = which( ! assays %in% assayNames(x))
		txt = paste("Assays are not found in dataset:", paste(head(assays[idx]), collapse=', '))
		stop(txt)
	}

	# extract metadata shared across assays
	data_constant = as.data.frame(data)

	pkeys = x@pkeys

	# for each assay
	resList = lapply( assays, function( k ){

		if( !quiet ) message('  ', k,'...', appendLF=FALSE)
		startTime = Sys.time()

		geneExpr = assay(x, k)

		# get names of samples to extract from metadata
		ids = colnames(geneExpr)

		# merge data_constant and pmetadata based on pkeys and assay k
		data2 = merge_metadata(data_constant[ids,,drop=FALSE], metadata(x), pkeys, k)
		data2 = droplevels(data2)

		# drop any constant terms from the formula
		form_mod = removeConstantTerms(formula, data2)

		# Drop variables in a redundant pair
		form_mod = dropRedundantTerms(form_mod, data2)

		# check if formula contains variables
		if( length(all.vars(form_mod)) > 0 ){
			# fit linear mixed model for each gene
			# TODO add , L=L
			res = fitExtractVarPartModel( geneExpr, form_mod, data2, BPPARAM=BPPARAM,...,quiet=TRUE)
		}else{
			res = data.frame()
		}

		if( !quiet ) message(format(Sys.time() - startTime, digits=2))

		list(df = res, formula = form_mod, n_retain = ncol(geneExpr))
	})
	# name each result by the assay name
	names(resList) = assays

	if( !quiet ) message("\n")

	# Convert results to DataFrame in vpDF
	vplst = lapply( names(resList), function(id){

		# get variation partition results
		df = resList[[id]]$df

		if( nrow(df) > 0){
			res = data.frame(assay = id, 
				gene = rownames(df), 
				data.frame(df))
		}else{
			res = data.frame()
		}
		res
	})
	names(vplst) = names(resList) 

	# Use smartbind in case a variable is droped from the analysis
	df = do.call(smartbind, vplst)
	if( nrow(df) > 0){
		df$assay = factor(df$assay, names(resList))
	}

	# extract details
	df_details = lapply( names(resList), function(id){

		data.frame( assay = id,
			n_retain = resList[[id]]$n_retain,
			formula = paste(as.character(resList[[id]]$formula), collapse=''),
			formDropsTerms = ! equalFormulas( resList[[id]]$formula, formula)	)
	})
	df_details = do.call(rbind, df_details)

	ndrop = sum(df_details$formDropsTerms)

	if( ndrop > 0){
		warning("Terms dropped from formulas for ", ndrop, " assays.\n Run details() on result for more information")
	}

	new("vpDF", DataFrame(df), df_details=df_details)
})











