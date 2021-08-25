# Gabriel Hoffman
# April 1, 2021
#
# dreamlet uses linear mixed models in dream to perform differential expression in single cell data
  
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

#' Subset with brackets
#'
#' Subset with brackets
#'
#' @param x dreamletProcessedData object
#' @param i indeces to extract
#'
#' @rdname extract-methods
#' @aliases [,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("[", signature(x="dreamletProcessedData"),
	function(x, i){   
		res = new("dreamletProcessedData", x@.Data[i], data = x@data)
		names(res) = names(x)[i]
		res
	}
)


setGeneric('assayNames', SummarizedExperiment::assayNames)
setGeneric('assay', SummarizedExperiment::assay)
setGeneric('colData', SummarizedExperiment::colData)

#' Get assayNames
#' 
#' Get assayNames
#' 
#' @param x dreamletProcessedData object
#' @param ... other arguments
#'
#' @rdname assayNames-methods
#' @aliases assayNames,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("assayNames", signature(x="dreamletProcessedData"),
	function(x, ...){   
		names(x)
	}
)
#' Get assay
#' 
#' Get assay
#' 
#' @param x dreamletProcessedData object
#' @param i number indicating index, or string indicating assay
#' @param withDimnames not used
#' @param ... other arguments
#'
#' @rdname assay-methods
#' @aliases assay,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("assay", signature(x="dreamletProcessedData"),
	function(x, i, withDimnames=TRUE,...){   
		x[[i]]
	}
)


#' Extract colData from dreamletProcessedData
#' 
#' Extract colData from dreamletProcessedData
#'
#' @param x A dreamletProcessedData object
#' @param ... other arguments
#' @export
setMethod("colData", "dreamletProcessedData",
	function(x,...){
		x@data
})



#' Show object
#' 
#' Show object
#' 
#' @param object dreamletProcessedData object
#'
#' @rdname show-methods
#' @aliases show,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("show", "dreamletProcessedData",
	function(object){
		print(object)
	}
)

#' Print object
#' 
#' Print object
#' 
#' @param x dreamletProcessedData object
#' @param ... other arguments
#' 
#' @importFrom utils head tail
#' @export
#' @rdname print-methods
#' @aliases print,dreamletProcessedData,dreamletProcessedData-method
setMethod("print", "dreamletProcessedData",
	function(x,...){

		cat('class: dreamletProcessedData\n')
		cat('assays(', length(x), '): ', sep='')

		if( length(x) < 4 ){
			txt = paste(names(x), collapse=' ')
		}else{
			txt = paste(paste(head(names(x), 2), collapse=' '),
				'...',
				paste(tail(names(x), 2), collapse=' '))
		}
		cat(txt)

		df_count = lapply(x, function(obj) dim(obj$geneExpr))
		df_count = do.call(rbind, df_count)

		cat('\nSamples\n min:', min(df_count[,2]), '\n max:', max(df_count[,2]))
		cat('\nGenes\n min:', min(df_count[,1]), '\n max:', max(df_count[,1]), '\n\n')
	}
)

# setGene
#' Extract a subset of samples
#'
#' Extract a subset of samples
#' 
#' @param x dreamletProcessedData
#' @param ids column names to retain
#'
#' @export
subsetSamples = function(x, ids){

	stopifnot( is(x, 'dreamletProcessedData'))

	# for each assay
	for(i in seq_len(length(x)) ){

		# intersect ids with column names
		include = intersect(ids, colnames(x[[i]]$geneExpr))

		# extract samples with these column names
		x[[i]]$geneExpr = x[[i]]$geneExpr[,include]
	}

	x
}





# extract table of cell counts from 'int_colData'
# of pseudobulks as returned by 'aggregateData'
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment int_colData
.n_cells <- function(x) {
    y <- int_colData(x)$n_cells
    if (is.null(y)) return(NULL)
    if (length(metadata(x)$agg_pars$by) == 2)
        y <- as.matrix(data.frame(y, check.names = FALSE))
    return(as.table(y))
}





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
#' @param min.count min.count used by \code{edgeR::filterByExpr}
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
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#'
#' @export
processOneAssay = function( y, formula, data, n.cells, min.cells = 10, isCounts = TRUE, normalize.method = 'TMM', min.count = 10, BPPARAM = bpparam(),...){

	# nCells = extract from y

	# samples to include of they have enough observed cells
	include = (n.cells >= min.cells)

	# if no samples are retained
	if( sum(include) == 0){
		return(NULL)
	}

	# subset expression and data
	y = y[,include,drop=FALSE] 
	data = droplevels(data[include,,drop=FALSE])

	# per sample weights based on cell counts in sceObj
	w = n.cells[include] #weights by cell types

	# convert vector of sample weights to full matrix
	# each gene is weighted the same
	weights = asMatrixWeights(w, dim=c(nrow(y), ncol(y)))

	if( isCounts ){

		# Get count data and normalize
    	y = suppressMessages(DGEList(y, remove.zeros = TRUE))
    	y = calcNormFactors(y, method = normalize.method )

    	# drop any constant terms from the formula
		formula = removeConstantTerms(formula, data)

		# get samples with enough cells
		# filter genes
		design = model.matrix( subbars(formula), data)
		keep = filterByExpr(y, design, min.count=min.count)
		# keep[] = TRUE

		# create EList object storing gene expression and sample weights
		obj = new("EList", list(	E 	= y[keep,],
		 						weights = weights[keep,,drop=FALSE]))

		# since the sample weights are already in y, don't need to 
		# explicitly consider them here.
		geneExpr = voomWithDreamWeights( obj, formula, data, BPPARAM=BPPARAM,..., save.plot=TRUE, quiet=TRUE)

		result = list(	geneExpr = geneExpr, 
						isCounts = isCounts) 
	}else{
	 	
		# only include genes that show variation,
		# and have at least 5 nonzero values
	 	include = apply(y, 1, function(x)
	 		(var(x) > 0) & (sum(x!=0) > 5)
	 		) 

		# if data is already log2 CPM
		# create EList object storing gene expression and sample weights
		geneExpr = new("EList", list(E=y[include,,drop=FALSE], weights = weights[include,,drop=FALSE]))


		result = list(	geneExpr = geneExpr, 
						isCounts = isCounts)
	}

	result 
}


# since precision weights are not used, use the trend in the eBayes step
# trend = TRUE



#' Processing SingleCellExperiment to dreamletProcessedData
#'
#' For raw counts, estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param sceObj SingleCellExperiment object 
#' @param formula regression formula for differential expression analysis
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param min.count min.count used by \code{edgeR::filterByExpr}
#' @param pmetadata sample-specific metadata the varies across cell types.  This is merged with \code{colData(sceObj)} for each assay to make variables accessable to the formula
#' @param pkeys array of two strings indicating sample identifier and cell type identifier columns in pmetadata
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel  
#' @importFrom SummarizedExperiment colData assays assay assayNames
#' @importFrom S4Vectors metadata as.data.frame
#'
#' @export
processAssays = function( sceObj, formula, min.cells = 10, isCounts=TRUE, normalize.method = 'TMM', min.count = 10, pmetadata=NULL, pkeys=NULL, BPPARAM = bpparam(),...){

	# checks
	stopifnot( is(sceObj, 'SingleCellExperiment'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data_constant = as.data.frame(colData(sceObj))

	# pseudo-metadata about each sample for each cell type
	# subset this metadata for use be each assay (i.e. cell type)
	use_pmeta = FALSE
	if( !is.null(pmetadata) ){
		# cast to data.frame
		pmetadata = as.data.frame(pmetadata)
		use_pmeta = TRUE
	}

	# check that pkeys 
	found = pkeys %in% colnames(pmetadata) 
	if( any(!found) ){
		stop("pkeys not found in pmetadata: ", paste(pkeys[!found], collapse=", "))
	}

	# for each assay
	resList = lapply( assayNames(sceObj), function(k){

		y = assay(sceObj, k)
		# n.cells = metadata(sceObj)$n_cells[k,colnames(y),drop=FALSE]
		n.cells = .n_cells(sceObj)[k,colnames(y),drop=FALSE]

		if( use_pmeta ){
			# merge data with pseudo-metadata for this cell type
			pmetadata_sub = pmetadata[pmetadata[[pkeys[2]]] == k, ]

			# mergeand make sure order is the same
			data = merge(data_constant, pmetadata_sub, by.x="row.names", by.y=pkeys[1])
			rownames(data) = data$Row.names
			data = data[rownames(data_constant),]

			# identical(rownames(data), colnames(y))
		}else{
			data = data_constant
		}

		# processing counts with voom or log2 CPM
		processOneAssay(y, formula, data, n.cells, min.cells, isCounts, normalize.method, min.count = min.count, BPPARAM=BPPARAM,...)
	})
	names(resList) = assayNames(sceObj)

	new("dreamletProcessedData", resList, data = data_constant)
}



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
	function( x, formula, data = colData(x), L.list=NULL, min.cells = 10, isCounts=TRUE, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	standardGeneric("dreamlet")
})








#' @importFrom variancePartition getContrast dream
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @export
#' @rdname dreamlet
#' @aliases dreamlet,dreamletProcessedData-method
setMethod("dreamlet", "dreamletProcessedData",
	function( x, formula, data = colData(x), L.list=NULL, min.cells = 10, isCounts=TRUE, robust=FALSE, normalize.method = 'TMM', BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))

	# extract metadata shared across assays
	data = as.data.frame(data)

	# for each assay
	resList = lapply( x, function( procData ){

		# get names of samples to extract from metadata
		ids = colnames(procData$geneExpr)
		data_sub = droplevels(data[ids,,drop=FALSE])

		# drop any constant terms from the formula
		form_mod = removeConstantTerms(formula, data_sub)
			
		if( !is.null(form_mod) ){
			# construct contrasts based on design matrix for this datset
			if( ! is.null(L.list) ){
				L = lapply(L.list, function(coeffs){
					getContrast(procData$geneExpr, form_mod, data_sub, coeffs)
				})
				L = do.call(cbind, L)
			}else{
				L = NULL
			}

			fit = tryCatch( {
				# fit linear (mixed) model for each gene			
				dream( procData$geneExpr, form_mod, data_sub, L = L, BPPARAM=BPPARAM,..., quiet=TRUE)
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
				fit = eBayes(fit, robust=robust, trend=!procData$isCounts)
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

	resList
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
	function( x, formula, data = colData(x), BPPARAM = bpparam(),...){

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
	function( x, formula, data = colData(x), BPPARAM = bpparam(),...){

	# checks
	# stopifnot( is(x, 'dreamletProcessedData'))
	stopifnot( is(formula, 'formula'))
	
	# extract metadata shared across assays
	data = as.data.frame(data)

	# for each assay
	resList = lapply( x, function( procData ){

		# get names of samples to extract from metadata
		ids = colnames(procData$geneExpr)
		data_sub = droplevels(data[ids,,drop=FALSE])

		# drop any constant terms from the formula
		form_mod = removeConstantTerms(formula, data_sub)

		# fit linear mixed model for each gene
		# TODO add , L=L
		fitExtractVarPartModel( procData$geneExpr, form_mod, data_sub, BPPARAM=BPPARAM,...,quiet=TRUE)
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














