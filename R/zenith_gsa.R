# Gabriel Hoffman
# April 20, 2021
# 
# Apply zenith to dreamlet result


#' Perform zenith analysis
#' 
#' Perform zenith analysis
#'
#' @param fit results from \code{dreamlet}
#' @param coefs coefficients to test using \code{topTable(fit, coef)}
#' @param geneSets a list of gene sets
#' @param n_genes_min minumum number of genes in a geneset
#'  
#' @rdname zenith_gsa-methods
#' @export
setGeneric('zenith_gsa', function(fit, coefs, geneSets, n_genes_min = 10){
	standardGeneric("zenith_gsa")
	})



#' @importFrom limma ids2indices
#' @importFrom zenith zenith recodeToList
#' @importFrom stats p.adjust
#'
#' @rdname zenith_gsa-methods
#' @aliases zenith_gsa,dreamletResult,ANY,GeneSetCollection-method
#' @export
setMethod("zenith_gsa", signature(fit="dreamletResult", coefs="ANY", geneSets = "GeneSetCollection"),
	function(fit, coefs, geneSets, n_genes_min = 10){

	# convert GeneSetCollection to list
	geneSets.lst = recodeToList( geneSets )

	# for each assay
	df_zenith = lapply( assayNames(fit), function(k){

		# extract dream fit
		fit_local = assay(fit, k)

		# Map from Ensembl genes in geneSets_GO to 
		# from trimmed Ensembl names from RNA-seq data 
		index = ids2indices( geneSets.lst, rownames(fit_local))
		   
		# filter by size of gene set
		index = index[sapply(index, length) >= n_genes_min]

		# for each coefficient selected
		df_res = lapply( coefs, function(coef){
			# run zenith on dream fits
			df_res = zenith(fit_local, coef, index)
			
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





