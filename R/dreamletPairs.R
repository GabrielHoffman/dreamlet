# Gabriel Hoffman
# Nov 11, 2021


# @importFrom edgeR

# min.count = 5
# normalize.method = "TMM"
# method = c("random", "fixed", "none"),



#' Differential expression between pair of assays
#'
#' Perform differential expression between a pair of assays using linear (mixed) models
#'
#' @param pb pseudobulk data as SingleCellExperiment object
#' @param assay.pairs array of two entries specifying assays (i.e. cell clusters) to compare
#' @param method account for repeated measures from donors using a "random" effect, a "fixed" effect, or "none"
#' @param formula covariates to include in the analysis.
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param useCountsWeights use cell count weights
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param min.count minimum number of reads for a gene to be consider expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param robust logical, use eBayes method that is robust to outlier genes
# @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param quiet show messages
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
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
#' 	assay = "counts",    
#' 	cluster_id = 'cluster_id', 
#' 	sample_id = 'sample_id',
#' 	verbose=FALSE)
#' 
#' # compare first two assays (i.e. cell types)
#' ct.pairs = assayNames(pb)[1:2]
#' 
#' # run comparison
#' # use method = 'fixed' here since it is faster
#' res = dreamletPairs( pb, ct.pairs, method="fixed")
#' 
#' head(res)
#' 
#' @importFrom variancePartition dream eBayes topTable
#' @export
dreamletPairs = function( pb, assay.pairs, method = c("random", "fixed", "none"), formula = ~1, min.cells = 10, useCountsWeights = TRUE, isCounts=TRUE, normalize.method = 'TMM', min.count = 10, robust=FALSE, quiet=FALSE, BPPARAM = SerialParam(),...){

	method = match.arg(method)

	# check that assay.pairs has two entries
	if( length(assay.pairs) != 2){
		stop("assay.pairs must have 2 entries")
	}

	# check that assay.pairs is present in pseudobulk
	test = assay.pairs %in% assayNames(pb)
	if( any(!test) ){
		txt = paste0("Assays not found: ", paste0(assay.pairs[!test], collapse=', '))
		stop(txt)
	}

    # extract pseudobulk counts for each cell type
	countsMatrix1 = assay(pb, assay.pairs[1])
	colnames(countsMatrix1) = paste0(assay.pairs[1], '_', colnames(countsMatrix1))

	countsMatrix2 = assay(pb, assay.pairs[2])
	colnames(countsMatrix2) = paste0(assay.pairs[2], '_', colnames(countsMatrix2))

	countsMatrix = cbind(countsMatrix1, countsMatrix2)

	# create new data object
	# add variable for cell cluster membership	
	data = as.data.frame(rbind(colData(pb), colData(pb)))
	data$Sample = factor(rep(rownames(colData(pb)), 2))
	data$cellCluster = c(rep(assay.pairs[1], ncol(countsMatrix1)), rep(assay.pairs[2], ncol(countsMatrix2)) )
	data$cellCluster = factor(data$cellCluster, assay.pairs)
	rownames(data) = colnames(countsMatrix)

	# create formula to evalute voom and differential expression
	form = switch(method, 
			'random' = {update.formula( formula, ~ cellCluster + (1|Sample))},
			'fixed' = {update.formula( formula, ~ cellCluster + Sample)},
			'none' = {update.formula( formula, ~ cellCluster)})

	# extract cell counts for two specified assays
	n.cells = c(cellCounts(pb)[,assay.pairs] )

	# processing counts with voom or log2 CPM
	vobj = processOneAssay(countsMatrix, form, data, n.cells, min.cells, isCounts, normalize.method, min.count = min.count, useCountsWeights=useCountsWeights, BPPARAM=BPPARAM,...)

	# perform differential expression regression analysis
	idx = rownames(data) %in% colnames(vobj$E)
	fit = dream( vobj, form, data[idx,], BPPARAM=BPPARAM,..., quiet=TRUE)

	# borrow information across genes with the Empirical Bayes step
	fit = eBayes(fit, robust=robust, trend=!vobj$isCounts)

	# extract results
	key = paste0('cellCluster', assay.pairs[2])
	res = topTable(fit, coef=key, number=Inf)

	res
}

# debug(compareClusterPairs)

# source("/Users/gabrielhoffman/workspace/repos/dreamlet/R/compareClusterPairs.R")

# res = dreamletPairs( pb,  assay.pairs=assayNames(pb)[1:2] )
# res1 = dreamletPairs( pb, assay.pairs=assayNames(pb)[1:2], method = "fixed" )
# res2 = dreamletPairs( pb, assay.pairs=assayNames(pb)[1:2], method = "none" )



# df = merge(res1, res, by="row.names")

# plot(df$t.x, df$t.y)
# abline(0, 1, col="red")


















