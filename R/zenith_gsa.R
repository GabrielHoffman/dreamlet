# Gabriel Hoffman
# April 20, 2021
# 
# Apply zenith to dreamlet result

#' zenith gene set analysis on dreamlet results
#' 
#' zenith gene set analysis on dreamlet results
#'
#' @param res.dl results from \code{dreamlet}
#' @param coef coefficient to test using \code{topTable(fit, coef)}
#' @param geneSets a list of gene sets
#' @param n_genes_min minumum number of genes in a geneset
#'  
#' @importFrom limma ids2indices
#' @importFrom zenith zenith
#' @importFrom stats p.adjust
#'
#' @export
zenith_gsa = function(res.dl, coef, geneSets, n_genes_min = 10){

	# for each assay
	df_zenith = lapply( names(res.dl), function(assay){

		# extract dream fit
		fit = res.dl[[assay]]

		# Map from Ensembl genes in geneSets_GO to 
		# from trimmed Ensembl names from RNA-seq data 
		index = ids2indices( geneSets, rownames(fit))
		   
		# filter by size of gene set
		index = index[sapply(index, length) >= n_genes_min]

		# run zenith on dream fits
		df_res = zenith(fit, coef, index)

		data.frame(Assay = assay, Geneset = rownames(df_res), df_res)
		})
	names(df_zenith) = names(res.dl)
	df_zenith = do.call(rbind, df_zenith)
	rownames(df_zenith) = c()

	# Compute FDR using p-values from all assays
	df_zenith$FDR = p.adjust(df_zenith$PValue, "BH")

	df_zenith
}
