# Gabriel Hoffman
# April 1, 2021
#
# dreamlet uses linear mixed models in dream to perform differential expression in single cell data

# local definition so methods in this file have this class
setClass("dreamletProcessedData", contains="list", slots = c(data = 'data.frame', metadata='data.frame', pkeys="vector"))

#' Class dreamletResult
#'
#' Class \code{dreamletResult} stores results produced by \code{dreamlet()} to give a standard interface for downstream analysis
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
#' @return show data stored in object
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
#' @return print data stored in object
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

		if( is.null(df_count) ){
			cat("No assays retained\n")
		}else{
			cat('Genes:\n min:', min(df_count[,1]), '\n max:', max(df_count[,1]), '\n')
		}

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
#' @param obj A \code{dreamletResult} object
#'
#' @return array storing names of coefficients
#'
#' @examples
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
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data 
#' res.dl = dreamlet( res.proc, ~ group_id)
#' 
#' # show coefficients estimated for each cell type
#' coefNames(res.dl)
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

	unique(c(unlist(lapply(obj, function(x) colnames(coef(x))) )))
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


	res = new("dreamletResult", x@.Data[i], 
				df_details = details(x)[i,,drop=FALSE])
	names(res) = names(x)[i]
	res

	}
)





#' Convert list of regression fits to \code{dreamletResult}
#'
#' Convert list of regression fits to \code{dreamletResult} for downstream analysis
#'
#' @param fitList list of regression fit with \code{dream()}
#' @param df_details \code{data.frame} storing assay details
#' 
#' @details Useful for combining multiple runs of \code{dreamletCompareClusters()} into a single \code{dreamletResult} for downstream analysis
#' @examples
#' library(muscat)
#' library(SingleCellExperiment) 
#' 
#' data(example_sce)
#' 
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce, 
#'   assay = "counts",    
#'   cluster_id = 'cluster_id', 
#'   sample_id = 'sample_id',
#'  verbose=FALSE)
#' 
#' # first comparison
#' ct.pairs =  c("B cells", "CD14+ Monocytes")
#' fit = dreamletCompareClusters( pb, ct.pairs, method="fixed")
#' 
#' # second comparison
#' ct.pairs2 =  c("B cells", "CD8 T cells")
#' fit2 = dreamletCompareClusters( pb, ct.pairs2, method="fixed")
#' 
#' # Make a list storing each result with a meaningful name
#' fitList = list()
#' 
#' id = paste0('[', ct.pairs[1], ']_vs_[', ct.pairs[2], ']')
#' fitList[[id]] = fit
#' 
#' id = paste0('[', ct.pairs2[1], ']_vs_[', ct.pairs2[2], ']')
#' fitList[[id]] = fit2
#' 
#' # create a dreamletResult form this list
#' res.compare = as.dreamletResult( fitList )
#' 
#' library(zenith)
#' go.gs = get_GeneOntology("CC", to="SYMBOL")
#' 
#' # Run zenith on each result.
#' # The coef 'compare' is used for each model
#' res_zenith = zenith_gsa(res.compare, go.gs, coef = 'compare')
#' 
#' @importFrom methods new is
#' @export
as.dreamletResult = function(fitList, df_details=NULL){

	# check that names are not empty
	if( any(names(fitList) == '') ){
		stop("names(fitList) must not contain empty names")
	}

	# check that names are unique and not empty
	if( length(unique(names(fitList))) != length(fitList) ){
		stop("names(fitList) must be unique")
	}

	# check that entries are the result of a dream() fit
	if( any(!sapply(fitList, is, class2="MArrayLM")) ){
		stop("Each element of fitList must be a fit from dream()")
	}

	if( is.null(df_details) ){
		res = new("dreamletResult", fitList)
	}else{
		res = new("dreamletResult", fitList, df_details = df_details)	
	}

	res
}




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
#' @return \code{data.frame} storing hypothesis test for each gene and cell type
#'
#' @examples
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
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data 
#' res.dl = dreamlet( res.proc, ~ group_id)
#' 
#' # show coefficients estimated for each cell type
#' coefNames(res.dl)
#' 
#' # extract results using limma-style syntax
#' # combines all cell types together
#' # adj.P.Val gives study-wide FDR   
#' topTable(res.dl, coef="group_idstim", number=3)
#' 
#' @seealso \code{limma::topTable()}, \code{variancePartition::topTable()}
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
		
		if( any(!coef %in% coefNames(fit)) ){
			stop("coef must be in coefNames")
		}

		# Run topTable on each assay
		res = lapply( assayNames(fit), function(k){
			fit1 = assay(fit, k)

			if( is.null(genelist) ) genelist = rownames(fit1)

			# if coef is found, and at least one entry of genelist in is fit1
			good = all(coef %in% colnames(coef(fit1))) & any(rownames(fit1) %in% genelist)

			if( good ){
				tab = topTable(fit1, coef = coef, number = Inf, genelist = genelist, sort.by = "none", p.value=p.value, lfc=lfc, confint=confint)
				if( nrow(tab) > 0 ){
					tab = tab[!is.na(tab$ID),]

					# if doesn't have z.std, add it
					if( ! "z.std" %in% colnames(tab) ){
						tab$z.std = tab$t
					}

					res = data.frame(assay = k, tab)
				}else{
					res = NULL
				}
			}else{
				res = NULL
			}
			res
		})
		# combine across assays
		res = DataFrame(do.call(rbind, res))

		# remove rownames
		rownames(res) = c()

		if( nrow(res) == 0){
			stop("No results were found matching given criteria")
		}

		# apply multiple testing across *all* tests
		# subset based on number afterwards
		res$adj.P.Val = p.adjust( res$P.Value, adjust.method)

		opt = c('logFC', 'AveExpr', 'P', 't', 'B', 'none')
		if( ! sort.by %in% opt){
			stop("sort.by must be in: ", paste0(opt, collapse=', '))
		}

		# sorting
		ord <- switch(sort.by, logFC = order(abs(res$logFC), decreasing = TRUE), 
			AveExpr = order(res$AveExpr, decreasing = TRUE), 
			P = order(res$P.Value, decreasing = FALSE), 
			t = order(abs(t), decreasing = TRUE), 
			B = order(res$B, decreasing = TRUE), 
			none = seq_len(nrow(res)) )

		head(res[ord,], number)
	}
)



#' Test if coefficient is different from a specified value
#'
#' Test if coefficient is different from a specified value
#'
#' @param fit dreamletResult object
#' @param lfc a minimum log2-fold-change below which changes not considered scientifically meaningful
#' @param coef which coefficient to test
#' @param number number of genes to return
#' @param sort.by column to sort by
#'
#' @return \code{DataFrame} storing hypothesis test for each gene and cell type
#'
#' @examples
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
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data 
#' res.dl = dreamlet( res.proc, ~ group_id)
#' 
#' # show coefficients estimated for each cell type
#' coefNames(res.dl)
#' 
#' # extract results using limma-style syntax
#' # combines all cell types together
#' # adj.P.Val gives study-wide FDR
#' getTreat(res.dl, coef="group_idstim", number=3)
#' 
#' @seealso \code{limma::topTreat()}, \code{variancePartition::getTreat()}
#' @importFrom variancePartition getTreat
#' @rdname getTreat-methods
#' @aliases getTreat,dreamletResult,dreamletResult-method
#' @export
setMethod("getTreat", signature(fit="dreamletResult"),
	function(fit, lfc=log2(1.2), coef=NULL, number=10, sort.by = "p"){

		if( any(!coef %in% coefNames(fit)) ){
			stop("coef must be in coefNames")
		}

		adjust.method = "BH"

		# Run topTable on each assay
		res = lapply( assayNames(fit), function(k){
			fit1 = assay(fit, k)

			# if coef is not found 
			if( all(coef %in% colnames(coef(fit1))) ){
				tab = getTreat(fit1, lfc = lfc,
									coef = coef, 
									number = Inf)
				res = data.frame(assay = k, tab)
			}else{
				res = NULL
			}
			res
		})
		# combine across assays
		res = DataFrame(do.call(rbind, res))

		# remove rownames
		rownames(res) = c()

		# apply multiple testing across *all* tests
		# subset based on number afterwards
		res$adj.P.Val = p.adjust( res$P.Value, adjust.method)

		opt = c('logFC', 'AveExpr', 'P', 't', 'B', 'none')
		if( ! tolower(sort.by) %in% tolower(opt)){
			stop("sort.by must be in: ", paste0(opt, collapse=', '))
		}

		# sorting
		ord <- switch( tolower(sort.by), 
			logfc = order(abs(res$logFC), decreasing = TRUE), 
			aveexpr = order(res$AveExpr, decreasing = TRUE), 
			p = order(res$P.Value, decreasing = FALSE), 
			t = order(abs(t), decreasing = TRUE), 
			b = order(res$B, decreasing = TRUE), 
			none = seq_len(nrow(res)) )

		head(res[ord,], number)
	}
)



#' Differential expression for each assay
#'
#' Perform differential expression for each assay using linear (mixed) models
#'
#' @param x SingleCellExperiment or dreamletProcessedData object 
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(x)}
#' @param contrasts character vector specifying contrasts specifying linear combinations of fixed effects to test.  This is fed into \code{makeContrastsDream( formula, data, contrasts=contrasts)}
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
# @param isCounts logical, indicating if data is raw counts
#' @param robust logical, use eBayes method that is robust to outlier genes
# @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param quiet show messages
#' @param BPPARAM parameters for parallel evaluation
#' @param use.eBayes should \code{eBayes} be used on result? (defualt: TRUE)
#' @param ... other arguments passed to \code{dream}
#'
#' @details
#' Fit linear (mixed) model on each cell type separately.  For advanced use of contrasts see \code{variancePartition::makeContrastsDream()} and vignette \url{https://gabrielhoffman.github.io/variancePartition/articles/dream.html#advanced-hypothesis-testing-1}.  
#'
#' @return Object of class \code{dreamletResult} storing results for each cell type
#'
#' @examples
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
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data 
#' res.dl = dreamlet( res.proc, ~ group_id)
#' 
#' # show coefficients estimated for each cell type
#' coefNames(res.dl)
#' 
#' # extract results using limma-style syntax
#' # combines all cell types together
#' # adj.P.Val gives study-wide FDR
#' topTable(res.dl, coef="group_idstim", number=3)
#' 
#' @import BiocParallel  
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @seealso \code{variancePartition::dream()}, \code{variancePartition::makeContrastsDream()}
#' @export
setGeneric("dreamlet", 
	function( x, formula, data = colData(x), assays = assayNames(x), contrasts=NULL, min.cells = 10, robust=FALSE, quiet=FALSE, BPPARAM = SerialParam(), use.eBayes=TRUE,...){

	standardGeneric("dreamlet")
})





#' @importFrom variancePartition getContrast dream makeContrastsDream
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @import Rdpack
#' @export
#' @rdname dreamlet
#' @aliases dreamlet,dreamletProcessedData-method
setMethod("dreamlet", "dreamletProcessedData",
	function( x, formula, data = colData(x), assays = assayNames(x), contrasts=NULL, min.cells = 10, robust=FALSE, quiet=FALSE, BPPARAM = SerialParam(), use.eBayes=TRUE,...){

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

	# remove samples with missing covariate data
	idx = sapply(all.vars(formula), function(v) {
	        which(is.na(data_constant[[v]]))
    })
    idx = unique(unlist(idx))    

    if( length(idx) > 1){
	    data_constant = droplevels(data_constant[-idx,,drop=FALSE])
	}

	pkeys = x@pkeys

	# for each assay
	resList = lapply( assays, function( k ){

		if( !quiet ) message('  ', k,'...', appendLF=FALSE)
		startTime = Sys.time()

		geneExpr = assay(x, k)

		# get names of samples to extract from 
		# intersecting between geneExpr and metadata
		ids = intersect(colnames(geneExpr), rownames(data_constant)) 
		geneExpr = geneExpr[,ids,drop=FALSE]

		# merge data_constant and pmetadata based on pkeys and assay k
		data2 = merge_metadata(data_constant[ids,,drop=FALSE], metadata(x), pkeys, k)
		data2 = droplevels(data2)

		# drop any constant terms from the formula
		form_mod = removeConstantTerms(formula, data2)
		
		# Drop variables in a redundant pair
		form_mod = dropRedundantTerms(form_mod, data2)

		# drop any constant terms from the formula
		if( length(all.vars(form_mod)) > 0 ){

			# get contrasts customized for the formula for this cell type
			if( ! is.null(contrasts) ){
				L = makeContrastsDream( form_mod, data2, contrasts=contrasts, nullOnError=TRUE)
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

				if( use.eBayes ){

					# use counts directly if is an EList
					isCounts = ifelse(is(geneExpr, "EList"), TRUE, FALSE)

					# borrow information across genes with the empirical Bayes step
					fit = eBayes(fit, robust=robust, trend=isCounts)
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
	names(resList) = assays

	if( !quiet ) message("\n")

	# extract fit
	fitList = lapply(resList, function(x) x$fit)

	# only keep entries that are not NULL
	# NUll is returned when coef of interest is dropped
	fitList = fitList[!vapply(fitList, is.null, FUN.VALUE=logical(1))]

	# extract details
	df_details = lapply( names(resList), function(id){

		data.frame( assay = id,
			n_retain = resList[[id]]$n_retain,
			formula = Reduce(paste, deparse(resList[[id]]$formula)),
			formDropsTerms = ! equalFormulas( resList[[id]]$formula, formula)	)
	})
	df_details = do.call(rbind, df_details)

	ndrop = sum(df_details$formDropsTerms)

	if( ndrop > 0){
		warning("Terms dropped from formulas for ", ndrop, " assays.\n Run details() on result for more information")
	}

	new("dreamletResult", fitList, df_details = df_details)
})



























