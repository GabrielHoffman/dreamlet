

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
#' @param min.count minimum number of reads for a gene to be consider expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @import BiocParallel 
#' @import limma 
#' @importFrom variancePartition voomWithDreamWeights
#' @importFrom edgeR calcNormFactors filterByExpr DGEList 
#' @importFrom methods is new
#' @importFrom stats model.matrix var
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @importFrom lme4 subbars
#'
#' @export
processOneAssay = function( y, formula, data, n.cells, min.cells = 10, isCounts = TRUE, normalize.method = 'TMM', min.count = 10, BPPARAM = SerialParam(),...){

    checkFormula( formula, data)
    if( is.null(n.cells) ){
    	stop("n_cells must not be NULL")
    }

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
		# design = model.matrix( subbars(formula), data)
		# Design often includes batch and donor, which are very small
		# 	this causes too many genes to be retained 
		keep = suppressWarnings(filterByExpr(y, min.count=min.count))

		# create EList object storing gene expression and sample weights
		obj = new("EList", list(	E 	= y[keep,],
		 						weights = weights[keep,,drop=FALSE]))

		# since the sample weights are already in y, don't need to 
		# explicitly consider them here.
		geneExpr = voomWithDreamWeights( obj, formula, data, BPPARAM=BPPARAM,..., save.plot=TRUE, quiet=TRUE)

		# save formula used after dropping constant terms
		geneExpr$formula = formula
	}else{
	 	
		# assumse already converted to log2 CPM

		# only include genes that show variation,
		# and have at least 5 nonzero values
	 	include = apply(y, 1, function(x)
	 		(var(x) > 0) & (sum(x!=0) > 5)
	 		) 

		# if data is already log2 CPM
		# create EList object storing gene expression and sample weights
		geneExpr = new("EList", list(E=y[include,,drop=FALSE], weights = weights[include,,drop=FALSE]))
	}

	geneExpr$isCounts = isCounts

	geneExpr 
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
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
processAssays = function( sceObj, formula, min.cells = 10, isCounts=TRUE, normalize.method = 'TMM', min.count = 10, pmetadata=NULL, pkeys=NULL, BPPARAM = SerialParam(),...){

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

		# check that pkeys 
		found = pkeys %in% colnames(pmetadata) 
		if( any(!found) ){
			stop("pkeys not found in pmetadata: ", paste(pkeys[!found], collapse=", "))
		}
	}else{
		pmetadata = data.frame()
		pkeys = array()
	}

	# for each assay
	resList = lapply( assayNames(sceObj), function(k){

		y = assay(sceObj, k)

		# dreamlet style
		n.cells = .n_cells(sceObj)[k,colnames(y),drop=FALSE]

		if(is.null(n.cells)){
			# muscat style
			n.cells = metadata(sceObj)$n_cells[k,colnames(y),drop=FALSE]
		}

		# merge data_constant and pmetadata based on pkeys and assay k
		data = merge_metadata(data_constant, pmetadata, pkeys, k)

		# processing counts with voom or log2 CPM
		processOneAssay(y, formula, data, n.cells, min.cells, isCounts, normalize.method, min.count = min.count, BPPARAM=BPPARAM,...)
	})
	names(resList) = assayNames(sceObj)

	new("dreamletProcessedData", resList, data = data_constant, metadata = pmetadata, pkeys=pkeys)
}


merge_metadata = function(data, mdata, pkeys, value){

	# if no mdata, return data
	if( nrow(mdata) == 0 ){
		return( data )
	}

	# subset mdata for this assay
	mdata_sub = mdata[mdata[[pkeys[2]]] == value, ]

	# merge and make sure order is the same
	dataOut = merge(data, mdata_sub, by.x="row.names", by.y=pkeys[1], all.x=TRUE)
	rownames(dataOut) = dataOut$Row.names
	dataOut[rownames(data),-1,drop=FALSE]
}







