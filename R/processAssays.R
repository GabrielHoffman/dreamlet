

#' Processing expression data from assay
#'
#' For raw counts, estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param y matrix of counts or log2 CPM
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param n.cells array of cell count for each sample
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param useCountsWeights use cell count weights
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
processOneAssay = function( y, formula, data, n.cells, min.cells = 10, isCounts = TRUE, normalize.method = 'TMM', min.count = 10, useCountsWeights = TRUE, BPPARAM = SerialParam(),...){

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

	# sample-level weights based on cell counts
	w_cells = n.cells[include] 

	if( ! useCountsWeights ){
		w_cells[] = 1
	}

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

		# weights from w_cells are used in calculating residuals that 
		# voom uses to compute precision weights
		geneExpr = voomWithDreamWeights( y[keep,], formula, data, weights = w_cells, BPPARAM=BPPARAM,..., save.plot=TRUE, quiet=TRUE)

		# combine empirical weights from voomWithDreamWeights
		# with weighting by the number of cells
		# w = w_cells - mean(w_cells)
		# geneExpr = applyQualityWeights(geneExpr, w_cells)

		# save formula used after dropping constant terms
		geneExpr$formula = formula
	}else{
	 	
		# assumes already converted to log2 CPM

		# only include genes that show variation,
		# and have at least 5 nonzero values
	 	include = apply(y, 1, function(x)
	 		(var(x) > 0) & (sum(x!=0) > 5)
	 		) 

		# if data is already log2 CPM
		# create EList object storing gene expression and sample weights
		geneExpr = new("EList", list(E=y[include,,drop=FALSE], weights = w_cells))

		geneExpr$formula = ~ 0
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
#' @param useCountsWeights use cell count weights
#' @param quiet show messages
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @return Object of class \code{dreamletProcessedData} storing voom-style normalized expression data
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
#' @import BiocParallel  
#' @importFrom SummarizedExperiment colData assays assay assayNames
#' @importFrom S4Vectors metadata as.data.frame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
processAssays = function( sceObj, formula, min.cells = 10, isCounts=TRUE, normalize.method = 'TMM', min.count = 10, pmetadata=NULL, pkeys=NULL, useCountsWeights=TRUE, quiet=FALSE, BPPARAM = SerialParam(),...){

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

		if( !quiet ) message('  ', k,'...', appendLF=FALSE)
		startTime = Sys.time()

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
		res = processOneAssay(y, formula, data, n.cells, min.cells, isCounts, normalize.method, min.count = min.count, useCountsWeights=useCountsWeights, BPPARAM=BPPARAM,...)

		if( !quiet ) message(format(Sys.time() - startTime, digits=2))

		res
	})
	names(resList) = assayNames(sceObj)

	# remove empty assays
	resList = resList[!vapply(resList, is.null, FUN.VALUE=logical(1))]

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







