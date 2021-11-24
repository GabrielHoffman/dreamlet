# Gabriel Hoffman
# April 20, 2021
# 
# Apply zenith to dreamlet result


#' Perform zenith analysis
#' 
#' Perform zenith analysis
#'
#' @param fit results from \code{dreamlet}
#' @param geneSets a list of gene sets
#' @param coefs coefficients to test using \code{topTable(fit, coef=coefs[i])}
#' @param use.ranks do a rank-based test \code{TRUE} or a parametric test \code{FALSE}? default: FALSE
#' @param n_genes_min minumum number of genes in a geneset
#' @param inter.gene.cor if NA, estimate correlation from data.  Otherwise, use specified value
#' @param progressbar if TRUE, show progress bar
#' @param statistic which statistic to use
#' @param ... other arguments
#'  
#' @return \code{data.frame} of results for each gene set and cell type 
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
#' # Load Gene Ontology database 
#' # use gene 'SYMBOL', or 'ENSEMBL' id
#' # use get_MSigDB() to load MSigDB
#' library(zenith)
#' go.gs = get_GeneOntology("CC", to="SYMBOL")
#'    
#' # Run zenith gene set analysis on result of dreamlet
#' res_zenith = zenith_gsa(res.dl, go.gs, 'group_idstim' )
#' 
#' # for each cell type select 5 genesets with largest t-statistic
#' # and 1 geneset with the lowest
#' # Grey boxes indicate the gene set could not be evaluted because
#' #    to few genes were represented
#' plotZenithResults(res_zenith, 5, 1)
#' 
#' @rdname zenith_gsa-methods
#' @export
setGeneric('zenith_gsa', function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE,...){
	standardGeneric("zenith_gsa")
	})



#' @importFrom limma ids2indices
#' @importFrom zenith zenith recodeToList
#' @importFrom stats p.adjust
#'
#' @rdname zenith_gsa-methods
#' @aliases zenith_gsa,dreamletResult,GeneSetCollection,ANY-method
#' @export
setMethod("zenith_gsa", signature(fit="dreamletResult", geneSets = "GeneSetCollection", coefs="ANY"),
	function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE,...){

	# convert GeneSetCollection to list
	geneSets.lst = recodeToList( geneSets )

	# for each assay
	df_zenith = lapply( assayNames(fit), function(k){

		# extract dream fit
		fit_local = assay(fit, k)

		# Map from genes to gene sets
		index = ids2indices( geneSets.lst, rownames(fit_local))
		   
		# filter by size of gene set
		index = index[vapply(index, length, FUN.VALUE=numeric(1)) >= n_genes_min]

		# for each coefficient selected
		df_res = lapply( coefs, function(coef){
			# run zenith on dream fits
			df_res = zenith(fit_local, coef, index, use.ranks=use.ranks, inter.gene.cor=inter.gene.cor, progressbar=progressbar)
			
			data.frame(assay = k, coef = coef, Geneset = rownames(df_res), df_res)
		})
		do.call(rbind, df_res)
	})
	names(df_zenith) = names(fit)
	df_zenith = do.call(rbind, df_zenith)
	rownames(df_zenith) = c()

	# Compute FDR using p-values from all assays
	df_zenith$FDR = p.adjust(df_zenith$PValue, "BH")

	df_zenith
})




#' @importFrom limma ids2indices
#' @importFrom zenith zenith recodeToList
#' @importFrom stats p.adjust
#'
#' @rdname zenith_gsa-methods
#' @aliases zenith_gsa,MArrayLM,GeneSetCollection,ANY-method
#' @export
setMethod("zenith_gsa", signature(fit="MArrayLM", geneSets = "GeneSetCollection", coefs="ANY"),
	function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE,...){

	# convert GeneSetCollection to list
	geneSets.lst = recodeToList( geneSets )

	# Map from genes to gene sets
	index = ids2indices( geneSets.lst, rownames(fit))
	   
	# filter by size of gene set
	index = index[vapply(index, length, FUN.VALUE=numeric(1)) >= n_genes_min]

	# for each coefficient selected
	df_zenith = lapply( coefs, function(coef){
		# run zenith on dream fits
		df_res = zenith(fit, coef, index, use.ranks=use.ranks, inter.gene.cor=inter.gene.cor, progressbar=progressbar)
		
		data.frame(coef = coef, Geneset = rownames(df_res), df_res)
	})
	df_zenith = do.call(rbind, df_zenith)
	
	df_zenith
})


#'
#' @importFrom limma ids2indices cameraPR
#' @importFrom zenith zenith recodeToList
#' @importFrom stats p.adjust
#' @importFrom ashr get_pm get_lfsr get_psd
#'
#' @rdname zenith_gsa-methods
#' @aliases zenith_gsa,dreamlet_mash_result,GeneSetCollection,ANY-method
#' @export
setMethod("zenith_gsa", signature(fit="dreamlet_mash_result", geneSets = "GeneSetCollection", coefs="ANY"),
	function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE, statistic=c('tstatistic', 'abs(tstatistic)', 'logFC', 'abs(logFC)')){

	statistic = match.arg(statistic)

	statMat = switch( statistic, 
		"logFC" = get_pm(fit$model),
		"abs(logFC)" = abs(get_pm(fit$model)),
		"tstatistic" = get_pm(fit$model) / get_psd(fit$model),
		"abs(tstatistic)" = abs(get_pm(fit$model) / get_psd(fit$model)))

	if( is.null(statMat) ){
		stop("statistic argument was not valid value: ", statistic)
	}

	# convert GeneSetCollection to list
	geneSets.lst = recodeToList( geneSets )

	# for each cell type (i.e. column)
	df_zenith = lapply( colnames(statMat), function(key){

		stats = statMat[,key,drop=TRUE]

		# remove NA values
		stats = stats[!is.na(stats)]

		# Map from genes to gene sets
		index = ids2indices( geneSets.lst, names(stats))
		   
		# filter by size of gene set
		index = index[vapply(index, length, FUN.VALUE=numeric(1)) >= n_genes_min]

		# run zenith on dream fits
		df_res = cameraPR(stats, index, use.ranks=use.ranks, inter.gene.cor=inter.gene.cor)
		
		data.frame(assay = key, coef = fit$coef, Geneset = rownames(df_res), df_res)
	})
	df_zenith = do.call(rbind, df_zenith)
	rownames(df_zenith) = c()

	# Compute FDR using p-values from all assays
	df_zenith$FDR = p.adjust(df_zenith$PValue, "BH")

	df_zenith
})








