# Gabriel Hoffman
# April 20, 2021
# 
# Apply zenith to dreamlet result

#' zenith gene set analysis on dreamlet results
#' 
#' zenith gene set analysis on dreamlet results
#'
#' @param x results from \code{dreamlet}
#' @param coefs coefficients to test using \code{topTable(fit, coef)}
#' @param geneSets a list of gene sets
#' @param n_genes_min minumum number of genes in a geneset
#'  
#' @importFrom limma ids2indices
#' @importFrom zenith zenith
#' @importFrom stats p.adjust
#'
#' @export
zenith_gsa = function(x, coefs, geneSets, n_genes_min = 10){

	# for each assay
	df_zenith = lapply( names(x), function(assay){

		# extract dream fit
		fit = x[[assay]]

		# Map from Ensembl genes in geneSets_GO to 
		# from trimmed Ensembl names from RNA-seq data 
		index = ids2indices( geneSets, rownames(fit))
		   
		# filter by size of gene set
		index = index[sapply(index, length) >= n_genes_min]

		# for each coefficient selected
		df_res = lapply( coefs, function(coef){
			# run zenith on dream fits
			df_res = zenith(fit, coef, index)
			
			data.frame(Assay = assay, GeneSet = rownames(df_res), df_res)
		})
		do.call(rbind, df_res)
	})
	names(df_zenith) = names(x)
	df_zenith = do.call(rbind, df_zenith)
	rownames(df_zenith) = c()

	# Compute FDR using p-values from all assays
	df_zenith$FDR = p.adjust(df_zenith$PValue, "BH")

	df_zenith
}





