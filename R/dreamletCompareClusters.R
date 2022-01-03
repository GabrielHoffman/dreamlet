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
#' @param min.count minimum number of reads for a gene to be consider expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param useCountsWeights use cell count weights
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
#'			baseline= c("CD4 T cells", "CD8 T cells"))
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

	# convert entries to strings
	assay.lst = lapply(assay.lst, as.character)

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

	n.samples = length(unique(data$Sample))
	n.cellCluster = length(unique(data$cellCluster))

	# processing counts with voom or log2 CPM
	vobj = processOneAssay(countsMatrix, form, data, n.cells, 
		min.cells = min.cells, 
		min.count = min.count,
		min.samples = min.samples,
		isCounts = isCounts,
		normalize.method = normalize.method,  
		quiet = quiet,
		useCountsWeights = useCountsWeights, 
		BPPARAM = BPPARAM,...)

	if( is.null(vobj) ){
		stop("No samples passed the filters.\n  Consider looser cutoffs for min.cells, min.count, min.samples")
	}

	if( nrow(vobj) < 4 ){
		stop("Only ", nrow(vobj), " samples passed the filters.\n  Consider looser cutoffs for min.cells, min.count, min.samples")
	}	

	# since vobj contains a subset of cells, also subset the data
	idx = match(colnames(vobj), rownames(data))
	data = data[idx,]

	# describe samples dropped by filtering	
	n.samples2 = length(unique(data$Sample))
	n.cellCluster2 = length(unique(data$cellCluster))

	if( ! quiet ){
		cat("Initial filtering...\n")
		if( n.samples2 - n.samples != 0){
			cat("  Dropped", (n.samples - n.samples2), '/', n.samples, "samples\n")
		}
		if( n.cellCluster2 - n.cellCluster != 0){
			cat("  Dropped", (n.cellCluster - n.cellCluster2), '/', n.cellCluster, "cell clusters\n")
		}
	}

	# If paired analysis is requested, and only one example of a Sample is found
	if( method %in% c("random", "fixed") ){

		# drop Samples and cellCluster with only 1 example, 
		# 	continue until no more changes are made 
		data2 = data

		n_remaining = nrow(data2)

		for(i in 1:100){

			tab = table(data2$Sample) > 1
			keep = data2$Sample %in% names(tab)[tab]
			data2 = data2[keep,]
			tab = table(data2$cellCluster) > 1
			keep = data2$cellCluster %in% names(tab)[tab]
			data2 = data2[keep,]
			
			if( nrow(data2) == n_remaining) break
			n_remaining = nrow(data2)
		}
		if(i == 100 || n_remaining == 0) stop("No samples remain after filtering")

		# dropped samples
		browser()

		# retain only these samples
		idx = match(rownames(data2), rownames(data))
		data = droplevels(data[idx,])
		vobj = vobj[,idx]

		# keep only retained cellCluster
		assay.lst$test = assay.lst$test[assay.lst$test %in% data$cellCluster]
		assay.lst$baseline = assay.lst$baseline[assay.lst$baseline %in% data$cellCluster]

		# describe samples dropped by filtering	on pairs
		n.samples3 = length(unique(data$Sample))
		n.cellCluster3 = length(unique(data$cellCluster))

		if( ! quiet ){
			cat("Filtering for paired samples...\n")
			if( n.samples3 - n.samples2 != 0){
				cat("  Dropped", (n.samples2 - n.samples3), '/', n.samples2, "samples\n")
			}
			if( n.cellCluster3 - n.cellCluster2 != 0){
				cat("  Dropped", (n.cellCluster2 - n.cellCluster3), '/', n.cellCluster2, "cell clusters\n")
			}
		}

		if( min(sapply(assay.lst, length)) == 0 ){
			stop("Insufficient cellClusters retained after filtering")
		}
	}

	# specify contrasts
	test = paste0('(', paste0('`cellCluster', assay.lst$test, '`', collapse=' + '), ') / ', length(assay.lst$test))
	baseline = 	paste0('(', paste0('`cellCluster', assay.lst$baseline, '`', collapse=' + '), ') / ', length(assay.lst$baseline))

	L = makeContrastsDream( form, data, contrasts = c(compare = paste(test, '-', baseline)))

	# perform differential expression regression analysis
	fit = dream( vobj, form, data, L=L, BPPARAM=BPPARAM,..., quiet=TRUE)

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


















