# Gabriel Hoffman
# April 21, 2021
# 
# plotZenithResults

#' Heatmap of zenith results
#'
#' Heatmap of zenith results showing genesets that have the top and bottom t-statistics from each assay.
#'
#' @param df result \code{data.frame} from \link{zenith_gsa}
#' @param ntop number of gene sets with highest t-statistic to show
#' @param nbottom number of gene sets with lowest t-statistic to show
#' 
#' @return Heatmap showing enrichment for gene sets and cell types
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
#' res_zenith = zenith_gsa(res.dl, coef = 'group_idstim', go.gs)
#' 
#' # for each cell type select 3 genesets with largest t-statistic
#' # and 1 geneset with the lowest
#' # Grey boxes indicate the gene set could not be evaluted because
#' #    to few genes were represented
#' plotZenithResults(res_zenith, 3, 1)
#' 
#' @importFrom reshape2 dcast
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom stats hclust dist
#' @importFrom ashr get_pm get_lfsr get_psd
#' 
#' @export
plotZenithResults = function(df, ntop=5, nbottom=5){

	delta = se = PValue = NULL 
	
	if( all(c('delta', 'se') %in% colnames(df)) ){
		# construct t-statistic from effect size and standard error
		df$tstat = with(df, delta/se)
	}else{
		df$tstat = with(df, qnorm(PValue, lower.tail=FALSE))
		df$tstat = df$tstat * ifelse(df$Direction == "Up", 1, -1)
	}

	# if 'assay' is not found
	if( is.na(match("assay", colnames(df))) ){
		df$assay = df$coef
	}

	# for each assay, return top and bottom genesets
	gs = lapply( unique(df$assay), function(assay){

		lapply( unique(df$coef), function(coef){

			# extract zenith results for one assay
			df_sub = df[(df$assay == assay) &(df$coef == coef), ]

			# sort t-statistics
			tstat_sort = sort(df_sub$tstat)

			cutoff1 = ifelse(nbottom > 0, tstat_sort[nbottom], -Inf)
			cutoff2 = ifelse(ntop > 0, tail(tstat_sort, ntop)[1], Inf)

			# keep genesets with highest and lowest t-statistics
			idx = (df_sub$tstat <= cutoff1) | (df_sub$tstat >= cutoff2) 

			df_sub$Geneset[idx]
		})
	})
	gs = unique(unlist(gs))

	# create matrix from retained gene sets
	M = dcast(df[df$Geneset %in% gs,], assay + coef ~ Geneset, value.var = "tstat")
	annot = M[,seq(1,2)]
	M = as.matrix(M[,-seq(1,2)])
	rownames(M) = annot$assay

	# Perform clustering on data in M
	success = tryCatch({
		# hcl1 <- hclust(dist(M))
		hcl2 <- hclust(dist(t(M)))
		TRUE
		}, error = function(e) FALSE)

	# if original clustering fails, 
	# replace NA's with 0
	if( ! success ){
		M_zero = M
		i = which(is.na(M_zero))
		if( length(i) > 0) M_zero[i] = 0
		# hcl1 <- hclust(dist(M_zero))
		hcl2 <- hclust(dist(t(M_zero)))
	}

	# set breaks
	zmax = max(abs(M), na.rm=TRUE)
	at = seq(0, round(zmax), length.out=3)
	at = sort(unique(c(-at, at)))

	# set colors
	col_fun = colorRamp2(c(-zmax, 0, zmax), c("blue", "white", "red"))

	# create heatmap
	hm = Heatmap(t(M),
		        name = "t-statistic", #title of legend
		        # column_title = "assay", row_title = "Gene sets",
		        row_names_gp = gpar(fontsize = 8),
			    width = nrow(M), 
			    height = ncol(M),
			    cluster_rows = hcl2,
				# cluster_columns = hcl1,
			    column_split = annot$coef,
			    cluster_column_slices = FALSE,
			    heatmap_legend_param = list(at = at,  direction = "horizontal", title_position="topcenter"), 
			    col = col_fun)

	draw(hm, heatmap_legend_side = "bottom")
}










