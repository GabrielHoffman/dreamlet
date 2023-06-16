

#' Processing expression data from assay
#'
#' For raw counts, filter genes and samples, then estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param y matrix of counts or log2 CPM
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param n.cells array of cell count for each sample
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count minimum number of reads for a gene to be considered expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param useCountsWeights use cell count weights
#' @param quiet show messages
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @return \code{EList} object storing log2 CPM and precision weights
#'
#' @seealso \code{processAssays()}
#' @importFrom BiocParallel SerialParam
#' @importClassesFrom limma EList
#' @importFrom variancePartition voomWithDreamWeights
#' @importFrom edgeR calcNormFactors filterByExpr DGEList 
#' @importFrom methods is new
#' @importFrom stats model.matrix var
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @importFrom lme4 subbars
#'
processOneAssay = function( y, formula, data, n.cells, min.cells = 10, min.count = 10, min.samples=4, min.prop = .4, isCounts = TRUE, normalize.method = 'TMM', useCountsWeights = TRUE, quiet = TRUE, BPPARAM = SerialParam(),...){

	checkFormula( formula, data)
	if( is.null(n.cells) ){
		stop("n_cells must not be NULL")
	}
	if( !is.matrix(y) ){
		y = as.matrix(y)
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

	# if there are too few remaining samples
	if( nrow(data) < min.samples ){
		return( NULL )
	}

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

		# Drop variables in a redundant pair
		formula = dropRedundantTerms(formula, data)

		# get samples with enough cells
		# filter genes
		# design = model.matrix( subbars(formula), data)
		# Design often includes batch and donor, which are very small
		# 	this causes too many genes to be retained 
		keep = suppressWarnings(filterByExpr(y, min.count=min.count, min.prop=min.prop))

		# weights from w_cells are used in calculating residuals that 
		# voom uses to compute precision weights
		# If model fails return NULL
		geneExpr = tryCatch(
			voomWithDreamWeights( y[keep,], formula, data, weights = w_cells, BPPARAM=BPPARAM,..., save.plot=TRUE, quiet=quiet),
			error = function(e) NULL)

		# save formula used after dropping constant terms
		if( !is.null(geneExpr) ) geneExpr$formula = formula

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

	if( !is.null(geneExpr) ) geneExpr$isCounts = isCounts

	geneExpr 
}


# since precision weights are not used, use the trend in the eBayes step
# trend = TRUE



#' Processing SingleCellExperiment to dreamletProcessedData
#'
#' For raw counts, estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells.  
#'
#' @param sceObj SingleCellExperiment object 
#' @param formula regression formula for differential expression analysis
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(sceObj)}
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count minimum number of reads for a gene to be considered expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts for a gene to be retained
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param useCountsWeights use cell count weights
#' @param quiet show messages
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @return Object of class \code{dreamletProcessedData} storing voom-style normalized expression data
#'
#' @details  For each cell cluster, samples with at least \code{min.cells} are retained. Only clusters with at least \code{min.samples} retained samples are kept. Genes are retained if they have at least \code{min.count} reads in at least \code{min.prop} fraction of the samples.  Current values are reasonable defaults, since genes that don't pass these cutoffs are very underpowered for differential expression analysis and only increase the multiple testing burden.  But values of \code{min.cells = 5} and \code{min.count = 5} are also reasonable if you want to include more genes in the analysis.
#'
#' The precision weights are estimated using the residuals fit from the specified formula.  These weights are robust to changes in the formula as long as the major variables explaining the highest fraction of the variance are included.
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
#' @importFrom BiocParallel SerialParam 
#' @importFrom S4Vectors metadata as.data.frame
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays assay
#'
#' @export
processAssays = function( sceObj, formula, assays = assayNames(sceObj), min.cells = 10, min.count = 10, min.samples=4, min.prop = .4, isCounts = TRUE, normalize.method = 'TMM', useCountsWeights = TRUE, quiet=FALSE, BPPARAM = SerialParam(),...){

	# checks
	stopifnot( is(sceObj, 'SingleCellExperiment'))
	stopifnot( is(formula, 'formula'))

	# check colnames of SCE
	if( is.null(colnames(sceObj)) ){
		stop("colnames(sceObj) is NULL.  Column names are needed for internal filtering")
	}

	# extract metadata shared across assays
	data_constant = droplevels(as.data.frame(colData(sceObj)))

	# check if assays are valid
	if( any( ! assays %in% assayNames(sceObj)) ){
		idx = which( ! assays %in% assayNames(sceObj))
		txt = paste("Assays are not found in dataset:", paste(head(assays[idx]), collapse=', '))
		stop(txt)
	}

	# extract cell counts
	n.cells_full = cellCounts(sceObj)

	# extract all unique colnames
	colNamesAll = unique(unlist(lapply(assayNames(sceObj), function(x) colnames(assay(sceObj, x)))))

	# check for colnames missing cell counts
	if( any(!colNamesAll %in% rownames(n.cells_full) ) ){
		stop("Cell counts could not be extracted.\n  Check that colnames(sceObj) or rownames(colData(sceObj))\n  have not been modified manually after running aggregateToPseudoBulk()")
	}	

	# for each assay
	resList = lapply( assays, function(k){

		if( !quiet ) message('  ', k,'...', appendLF=FALSE)
		startTime = Sys.time()

		y = assay(sceObj, k)

		# cell counts
		n.cells = n.cells_full[colnames(y),k,drop=FALSE]

		# merge data_constant (data constant for all cell types)
		# with metadata(sceObj)$aggr_means (data that varies)
		data = merge_metadata(data_constant, 
								metadata(sceObj)$aggr_means, 
								k,
								metadata(sceObj)$agg_pars$by)

		# processing counts with voom or log2 CPM
		res = processOneAssay(y, formula, data, n.cells,
				min.cells = min.cells,
				min.count = min.count,
				min.samples = min.samples, 
				min.prop = min.prop,
				isCounts = isCounts,
				normalize.method = normalize.method, 
				useCountsWeights = useCountsWeights, 
				BPPARAM = BPPARAM,...)

		if( !quiet ) message(format(Sys.time() - startTime, digits=2))

		res
	})
	names(resList) = assays

	exclude = vapply(resList, is.null, FUN.VALUE=logical(1))

	warning("Not enough samples retained or model fit fails: ", paste(names(resList)[exclude], collapse=", "))

	# remove empty assays
	resList = resList[!exclude]

	new("dreamletProcessedData", 
		resList, 
		data = data_constant, 
		metadata = metadata(sceObj)$aggr_means, 
		by = metadata(sceObj)$agg_pars$by)
}

# merge data_constant (data constant for all cell types)
# with metadata(sceObj)$aggr_means (data that varies)
#' @importFrom dplyr filter_at
merge_metadata = function(dataIn, md, cellType, by){ 

	# PASS R CMD check
	cell = NULL

	data = merge(dataIn, 
			dplyr::filter_at( md, by[1], ~ . == cellType),
			by.x = "row.names",
			by.y = by[2])
	rownames(data) = data$Row.names
	data = data[rownames(dataIn),,drop=FALSE]
	droplevels(data)
}










