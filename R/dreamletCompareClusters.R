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
#' @param assays array of two entries specifying assays (i.e. cell clusters) to compare, of a list of two sets of assays.
#' @param method account for repeated measures from donors using a "random" effect, a "fixed" effect, or "none"
#' @param formula covariates to include in the analysis.
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param useCountsWeights use cell count weights
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param useCountsWeights use cell count weights
#' @param min.count minimum number of reads for a gene to be consider expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param robust logical, use eBayes method that is robust to outlier genes
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
#' ct.pairs =  c("B cells", "CD14+ Monocytes")
#' 
#' # run comparison
#' # use method = 'fixed' here since it is faster
#' res = dreamletCompareClusters( pb, ct.pairs, method="fixed")
#' 
#' # genes with highest logFC are most highly expressed in 
#' # B cells compared to CD14+ Monocytes
#' head(res)
#'
#' # compare B cells versus the rest of the cell types
#' # 'rest' is a keyword indicator all other assays	
#' res = dreamletCompareClusters( pb, c("B cells", 'rest'), method="fixed")
#' 
#' # genes with highest logFC are most highly expressed in 
#' # B cells compared to all others
#' head(res)
#'
#' lst = list( test = c("CD14+ Monocytes", "FCGR3A+ Monocytes"), 
#'			baseline = c("CD4 T cells", "CD8 T cells"))
#'
#' # compare 2 monocyte clusters to two T cell clusters
#' res = dreamletCompareClusters( pb, lst, method="fixed")
#' 
#' # genes with highest logFC are most highly expressed in 
#' # monocytes compared to T cells
#' head(res)
#' 
#' @importFrom variancePartition dream eBayes topTable makeContrastsDream
#' @export
dreamletCompareClusters = function( pb, assays, method = c("random", "fixed", "none"), formula = ~1, min.cells = 10, min.count = 10, min.samples=4, isCounts = TRUE, normalize.method = 'TMM', useCountsWeights=TRUE, robust=FALSE, quiet=FALSE, BPPARAM = SerialParam(),...){

	method = match.arg(method)

	if( is.vector(assays) & ! is.list(assays) ){

		# check that assays has two entries
		if( length(assays) != 2){
			stop("assays must have 2 entries")
		}

		if( assays[1] == "rest" ){
			stop("assay 'rest' is only valid as the second entry")
		}

		# convert array to list with two entries
		if( 'rest' %in% assays ){
			j = which(assayNames(pb) == assays[1])
			assay.lst = list(test = assays[1], baseline=assayNames(pb)[-j])
		}else{
			assay.lst = list(test = assays[1], baseline=assays[2])
		}
	}else if( is.list(assays) ){
		assay.lst = assays
	}else{
		stop("assays argument must be vector or list")
	}

	# check that assay.lst has two entries
	if( length(assay.lst) != 2){
		stop("assay.lst must have 2 entries")
	}

	# check that all are valid assays
	if( any(!unlist(assay.lst) %in% assayNames(pb)) ){
		
		i = unlist(assay.lst) %in% assayNames(pb)

		txt = paste(unlist(assay.lst)[!i], collapse=', ')
		stop("Specified assays are not valid: ", txt)
	}

	# check that there are no overlapping assays
	shared = intersect(assay.lst[[1]], assay.lst[[2]])
	if( length(shared) > 0 ){

		txt = paste(shared, collapse=', ')

		stop("Specified assays shared between two groups: ", txt)
	}

    # extract merge pseudobulk counts for specified assays
	res = lapply( unlist(assay.lst), function(clstr){
		geneCounts = assay(pb, clstr)
		colnames(geneCounts) = paste0(clstr[1], '_', colnames(geneCounts))
		
		df = as.data.frame(colData(pb))
		df$cellCluster = clstr
		df$Sample = rownames(df)
		rownames(df) = colnames(geneCounts)

		list(geneCounts = geneCounts, df = df)
		})

	countsMatrix = do.call(cbind, lapply(res, function(x) x$geneCounts))
	data = do.call(rbind, lapply(res, function(x) x$df))
	rownames(data) = colnames(countsMatrix)

	# create formula to evalute voom and differential expression
	form = switch(method, 
			'random' = {update.formula( formula, ~ 0 + cellCluster + (1|Sample))},
			'fixed' = {update.formula( formula, ~ 0 + cellCluster + Sample)},
			'none' = {update.formula( formula, ~ 0 + cellCluster)})

	# extract cell counts for two specified assays
	n.cells = c(cellCounts(pb)[,unlist(assay.lst)] )

	# processing counts with voom or log2 CPM
	vobj = processOneAssay(countsMatrix, form, data, n.cells, 
		min.cells = min.cells, 
		min.count = min.count,
		min.samples = min.samples,
		isCounts = isCounts,
		normalize.method = normalize.method,  
		robust = robust, 
		quiet = quiet,
		useCountsWeights = useCountsWeights, 
		BPPARAM = BPPARAM,...)

	# specify contrasts
	test = paste0('(', paste0('`cellCluster', assay.lst$test, '`', collapse=' + '), ') / ', length(assay.lst$test))
	baseline = 	paste0('(', paste0('`cellCluster', assay.lst$baseline, '`', collapse=' + '), ') / ', length(assay.lst$baseline))

	L = makeContrastsDream( form, data, contrasts = c(compare = paste(test, '-', baseline)))

	# perform differential expression regression analysis
	idx = match(colnames(vobj), rownames(data))
	fit = dream( vobj, form, data[idx,], L=L, BPPARAM=BPPARAM,..., quiet=TRUE)

	# borrow information across genes with the Empirical Bayes step
	fit = eBayes(fit, robust=robust, trend=!vobj$isCounts)

	# extract results
	res = topTable(fit, coef='compare', number=Inf)

	res[order(res$t, decreasing=TRUE),]
}

# debug(compareClusterPairs)

# source("/Users/gabrielhoffman/workspace/repos/dreamlet/R/compareClusterPairs.R")

# res = dreamletPairs( pb,  assays=assayNames(pb)[1:2] )
# res1 = dreamletPairs( pb, assays=assayNames(pb)[1:2], method = "fixed" )
# res2 = dreamletPairs( pb, assays=assayNames(pb)[1:2], method = "none" )



# df = merge(res1, res, by="row.names")

# plot(df$t.x, df$t.y)
# abline(0, 1, col="red")


















