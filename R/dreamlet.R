# Gabriel Hoffman
# April 1, 2021
#
# dreamlet uses linear mixed models in dream to perform differential expression in single cell data

# local definition so methods in this file have this class
setClass("dreamletProcessedData", contains="list", slots = c(data = 'data.frame', metadata='data.frame', pkeys="vector"))

#' Class dreamletResult
#'
#' Class \code{dreamletResult} 
#'
#' @name dreamletResult-class
#' @rdname dreamletResult-class
#' @exportClass dreamletResult
setClass("dreamletResult", contains="list", slots=c(df_details = "data.frame"))



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
	        nms <- character(length(assays(x, withDimnames=FALSE)))
	    coolcat("assays(%d): %s\n", nms)

		df_count = lapply(x, function(obj) nrow(obj$coefficients))
		df_count = do.call(rbind, df_count)

		cat('Genes:\n min:', min(df_count[,1]), '\n max:', max(df_count[,1]), '\n')

		# metadata
	    nms <- names(details(x))
	    if (is.null(nms))
	        nms <- character(length(metadata(x, withDimnames=FALSE)))
	    coolcat("details(%d): %s\n", nms)

	    # show coef names
	    coolcat("coefNames(%d): %s\n", coefNames(x))
	}
)



#' Get coefficient names
#' 
#' Get coefficient names
#'
#' @param obj A dreamletResult object
#'
#' @rdname coefNames-methods
#' @export
setGeneric('coefNames', function(obj){
	standardGeneric("coefNames")
	})

#' @export
#' @rdname coefNames-methods
#' @aliases coefNames,dreamletResult-method
#' @importFrom stats coef
setMethod("coefNames", "dreamletResult",
	function(obj){		

	unique(c(unlist(sapply(obj, function(x) colnames(coef(x))))))
})






setGeneric('assayNames', SummarizedExperiment::assayNames)
setGeneric('assay', SummarizedExperiment::assay)
# setGeneric('colData', SummarizedExperiment::colData)
# setGeneric('metadata', S4Vectors::metadata)

#' Get assayNames
#' 
#' Get assayNames
#' 
#' @param x dreamletResult object
#' @param ... other arguments
#'
#' @rdname assayNames-methods
#' @aliases assayNames,dreamletResult,dreamletResult-method
#' @export
setMethod("assayNames", signature(x="dreamletResult"),
	function(x, ...){   
		names(x)
	}
)

#' Get assay
#' 
#' Get assay
#' 
#' @param x dreamletResult object
#' @param i number indicating index, or string indicating assay
#' @param withDimnames not used
#' @param ... other arguments
#'
#' @rdname assay-methods
#' @aliases assay,dreamletResult,dreamletResult-method
#' @export
setMethod("assay", signature(x="dreamletResult"),
	function(x, i, withDimnames=TRUE,...){   
		x[[i]]
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
#' @param coef coef
#' @param number number
#' @param genelist genelist
#' @param adjust.method adjust.method
#' @param sort.by sort.by
#' @param resort.by resort.by
#' @param p.value p.value
#' @param lfc lfc
#' @param confint confint
#'
#' @rdname topTable-methods
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
		res = lapply( assayNames(fit), function(k){
			fit1 = assay(fit, k)

			if( is.null(genelist) ) genelist = rownames(fit1)

			tab = topTable(fit1, coef = coef, number = Inf, genelist = genelist, p.value=p.value, lfc=lfc, confint=confint)
			data.frame(assay = k, tab)
		})
		# combine across assays
		res = DataFrame(do.call(rbind, res))

		# remove rownames
		rownames(res) = c()

		# apply multiple testing across *all* tests
		# subset based on number afterwards
		res$adj.P.Val = p.adjust( res$P.Value, adjust.method)

		# sorting
		ord <- switch(sort.by, logFC = order(abs(res$logFC), decreasing = TRUE), 
			AveExpr = order(res$AveExpr, decreasing = TRUE), 
			P = order(res$P.Value, decreasing = FALSE), 
			t = order(abs(t), decreasing = TRUE), 
			B = order(res$B, decreasing = TRUE), none = seq_len(nrow(res)) )

		head(res[ord,], number)
	}
)






#' Differential expression for each assay
#'
#' Perform differential expression for each assay using linear mixed models
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param contrasts character vector specifying contrasts specifying linear combinations of fixed effects to test
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param isCounts logical, indicating if data is raw counts
#' @param robust logical, use eBayes method that is robust to outlier genes
# @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param quiet show messages
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel  
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @seealso \code{dream}, \code{makeContrastsDream}
#'
#' @export
setGeneric("dreamlet", 
	function( x, formula, data = colData(x), contrasts=NULL, min.cells = 10, isCounts=TRUE, robust=FALSE, quiet=FALSE, BPPARAM = SerialParam(),...){

	standardGeneric("dreamlet")
})




#' @importFrom variancePartition getContrast dream
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @export
#' @rdname dreamlet
#' @aliases dreamlet,dreamletProcessedData-method
setMethod("dreamlet", "dreamletProcessedData",
	function( x, formula, data = colData(x), contrasts=NULL, min.cells = 10, isCounts=TRUE, robust=FALSE, quiet=FALSE, BPPARAM = SerialParam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data_constant = as.data.frame(data)

	pkeys = x@pkeys

	# for each assay
	resList = lapply( names(x), function( k ){

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

		# drop any constant terms from the formula
		if( length(all.vars(form_mod)) > 0 ){

			# get contrasts customized for the formula for this cell type
			if( ! is.null(contrasts) ){
				L = makeContrastsDream( form_mod, data2, contrasts=contrasts)
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

				if( ! missing(omiteBayes) ){
					# borrow information across genes with the Empircal Bayes step
					fit = eBayes(fit, robust=robust, trend=!geneExpr$isCounts)
				}
			}else{	
				fit = NULL
			}
		}else{
			fit = NULL
		}

		if( !quiet ) message(format(Sys.time() - startTime, digits=2))

		list(fit = fit, formula = form_mod, n_retain = ncol(geneExpr))
	})
	# name each result by the assay name
	names(resList) = names(x)

	if( !quiet ) message("\n")

	# extract fit
	fitList = lapply(resList, function(x) x$fit)

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

	new("dreamletResult", fitList, df_details = df_details)
})



























