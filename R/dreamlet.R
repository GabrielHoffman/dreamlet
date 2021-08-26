# Gabriel Hoffman
# April 1, 2021
#
# dreamlet uses linear mixed models in dream to perform differential expression in single cell data

#
# need to adapt voomByDreamWeights to allow per sample weights
# see use of sample weights here: 
# https://rdrr.io/bioc/limma/src/R/voomWithQualityWeights.R




#' Class dreamletResult
#'
#' Class \code{dreamletResult} 
#'
#' @name dreamletResult-class
#' @rdname dreamletResult-class
#' @exportClass dreamletResult
setClass("dreamletResult", contains="list")



#' Show object
#' 
#' Show object
#' 
#' @param object dreamletResult object
#'
#' @rdname show-methods
#' @aliases show,dreamletResult,dreamletResult-method
#' @export
setMethod("show", "dreamletResult",
	function(object){
		print(object)
	}
)

#' Print object
#' 
#' Print object
#' 
#' @param x dreamletResult object
#' @param ... other arguments
#' 
#' @importFrom utils head tail
#' @importFrom S4Vectors coolcat
#' @export
#' @rdname print-methods
#' @aliases print,dreamletResult,dreamletResult-method
setMethod("print", "dreamletResult",
	function(x,...){

		cat('class:', class(x), '\n')

		# assay
	    nms <- names(x)
	    if (is.null(nms))
	        nms <- character(length(assays(object, withDimnames=FALSE)))
	    coolcat("assays(%d): %s\n", nms)

		df_count = lapply(x, function(obj) nrow(obj$coefficients))
		df_count = do.call(rbind, df_count)

		cat('Genes:\n min:', min(df_count[,1]), '\n max:', max(df_count[,1]), '\n')
	}
)



#' Subset with brackets
#'
#' Subset with brackets
#'
#' @param x dreamletResult object
#' @param i indeces to extract
#'
#' @rdname extract-methods
#' @aliases [,dreamletResult,dreamletResult-method
#' @export
setMethod("[", signature(x="dreamletResult"),
	function(x, i){   
		new("dreamletResult", x[i])
	}
)


#' Table of Top Genes from dreamlet fit
#'
#' Extract a table of the top-ranked genes from a dreamlet fit.
#'
#' @param fit dreamletResult object
#' @param ... arguments passed to topTable
#'
#' @rdname extract-methods
#' @aliases topTable,dreamletResult,dreamletResult-method
#' @export
setMethod("topTable", signature(fit="dreamletResult"),
	function(fit,       
		coef = NULL,
       number = 10,
       genelist = NULL,
       adjust.method = "BH",
       sort.by = "P",
       resort.by = NULL,
       p.value = 1,
       lfc = 0,
       confint = FALSE){   
		
		# Run topTable on each assay
		res = lapply( names(x), function(k){
			fit = x[[k]]

			if( is.null(genelist) ) genelist = rownames(fit)

			tab = topTable(fit, coef = coef, number = number, genelist = genelist, p.value=p.value, lfc=lfc, confint=confint)
			data.frame(assay = k, tab)
		})
		# combine across assays
		res = DataFrame(do.call(rbind, res))

		# apply multiple testing across all tests
		res$adj.P.Val = p.adjust( res$P.Value, adjust.method)

		# sorting
		ord <- switch(sort.by, logFC = order(abs(res$logFC), decreasing = TRUE), 
			AveExpr = order(res$AveExpr, decreasing = TRUE), 
			P = order(res$P.Value, decreasing = FALSE), 
			t = order(abs(t), decreasing = TRUE), 
			B = order(res$B, decreasing = TRUE), none = seq_len(nrow(res)) )

		res[ord,] 
	}
)






#' Differential expression for each assay
#'
#' Perform differential expression for each assay using linear mixed models
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param L.list list of contrasts specifying linear combinations of fixed effects to tests
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param isCounts logical, indicating if data is raw counts
#' @param robust logical, use eBayes method that is robust to outlier genes
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel  
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#'
#' @export
setGeneric("dreamlet", 
	function( x, formula, data = colData(x), L.list=NULL, min.cells = 10, isCounts=TRUE, robust=FALSE, normalize.method = 'TMM', BPPARAM = SerialParam(),...){

	standardGeneric("dreamlet")
})



#' @importFrom variancePartition getContrast dream
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @export
#' @rdname dreamlet
#' @aliases dreamlet,dreamletProcessedData-method
setMethod("dreamlet", "dreamletProcessedData",
	function( x, formula, data = colData(x), L.list=NULL, min.cells = 10, isCounts=TRUE, robust=FALSE, normalize.method = 'TMM', BPPARAM = SerialParam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data_constant = as.data.frame(data)

	pkeys = x@pkeys

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

		if( !is.null(form_mod) ){
			# construct contrasts based on design matrix for this datset
			if( ! is.null(L.list) ){
				L = lapply(L.list, function(coeffs){
					getContrast(geneExpr, form_mod, data2, coeffs)
				})
				L = do.call(cbind, L)
			}else{
				L = NULL
			}

			fit = tryCatch( {
				# fit linear (mixed) model for each gene			
				dream( geneExpr, form_mod, data2, L = L, BPPARAM=BPPARAM,..., quiet=TRUE)
				}, 
				error = function(e) NULL)

			# if model is degenerate
			if( !is.null(fit) && ! any(is.na(fit$sigma)) ){

				if( !is.null(fit$rdf)){
					# keep genes with residual degrees of freedom > 1
					# this prevents failures later
					keep = which(fit$rdf >= 1)

					fit = fit[keep,]
				}

				# borrow information across genes with the Empircal Bayes step
				fit = eBayes(fit, robust=robust, trend=!geneExpr$isCounts)
			}else{	
				fit = NULL
			}
		}else{
			fit = NULL
		}
		fit
	})
	# name each result by the assay name
	names(resList) = names(x)

	new("dreamletResult", resList)
})



























