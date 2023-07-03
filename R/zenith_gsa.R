# Gabriel Hoffman
# April 20, 2021
# 
# Apply zenith to dreamlet result



setGeneric('zenith_gsa', zenith::zenith_gsa)

#' Perform gene set analysis using zenith
#' 
#' Perform a competitive gene set analysis accounting for correlation between genes.
#'
#' @param fit results from \code{dreamlet()}
#' @param geneSets \code{GeneSetCollection} 
#' @param coefs coefficients to test using \code{topTable(fit, coef=coefs[i])}
#' @param use.ranks do a rank-based test \code{TRUE} or a parametric test \code{FALSE}? default: FALSE
#' @param n_genes_min minimum number of genes in a geneset
#' @param inter.gene.cor if NA, estimate correlation from data.  Otherwise, use specified value
#' @param progressbar if TRUE, show progress bar
# @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments
#'  
#' @return \code{data.frame} of results for each gene set and cell type 
#'
#' @details This code adapts the widely used \code{camera()} analysis \insertCite{wu2012camera}{zenith} in the \code{limma} package \insertCite{ritchie2015limma}{zenith} to the case of linear (mixed) models used by \code{variancePartition::dream()}.
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
#' res_zenith = zenith_gsa(res.dl, go.gs, 'group_idstim', progressbar=FALSE )
#' 
#' # for each cell type select 3 genesets with largest t-statistic
#' # and 1 geneset with the lowest
#' # Grey boxes indicate the gene set could not be evaluted because
#' #    to few genes were represented
#' plotZenithResults(res_zenith, 3, 1)
#' 
#' @importFrom limma ids2indices
#' @importFrom zenith zenith 
#' @importFrom GSEABase geneIds
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply
#'
#' @rdname zenith_gsa-methods
#' @aliases zenith_gsa,dreamletResult,GeneSetCollection,ANY-method
#' @export
setMethod("zenith_gsa", signature(fit="dreamletResult", geneSets = "GeneSetCollection", coefs="ANY"),
	function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE,...){

	# check that coefs are in the dreamlet result
	if( any(!coefs %in% coefNames(fit)) ){
		i = which(!coefs %in% coefNames(fit))
		txt = paste0("coefs are not found in dreamletResult: ", paste(coefs[i], sep=','))
		stop(txt)
	}

	# convert GeneSetCollection to list
	geneSets.lst = geneIds( geneSets )
	rm(geneSets)

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
			# if coef is available
			if( (coef %in% colnames(coef(fit_local))) & (length(index) > 0) ){
				# run zenith on dream fits
				df_res = zenith(fit_local, coef, index, use.ranks=use.ranks, inter.gene.cor=inter.gene.cor, progressbar=progressbar)
			
				df = data.frame(assay = k, coef = coef, Geneset = rownames(df_res), df_res)
			}else{
				df = NULL
			}
			df
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







#' @return \code{data.frame} of results for each gene set and cell type 
#' @importFrom limma ids2indices cameraPR
#' @importFrom zenith zenith
#' @importFrom GSEABase geneIds
#' @importFrom stats p.adjust
#' @importFrom ashr get_pm get_lfsr get_psd
#'
#' @rdname zenith_gsa-methods
#' @aliases zenith_gsa,dreamlet_mash_result,GeneSetCollection,ANY-method
#' @export
setMethod("zenith_gsa", signature(fit="dreamlet_mash_result", geneSets = "GeneSetCollection", coefs="ANY"),
	function(fit, geneSets, coefs, use.ranks=FALSE, n_genes_min = 10, inter.gene.cor=0.01, progressbar=TRUE,...){
		
	args <- list(...)
	if( 'statistic' %in% names(args) ){
		statistic = args$statistic
	}else{
		statistic = "tstatistic"
	}

	statMat = switch( statistic, 
		"logFC" = get_pm(fit$model),
		"abs(logFC)" = abs(get_pm(fit$model)),
		"tstatistic" = get_pm(fit$model) / get_psd(fit$model),
		"abs(tstatistic)" = abs(get_pm(fit$model) / get_psd(fit$model)))

	if( is.null(statMat) ){
		stop("statistic argument was not valid value: ", statistic)
	}

	# convert GeneSetCollection to list
	geneSets.lst = geneIds( geneSets )

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








