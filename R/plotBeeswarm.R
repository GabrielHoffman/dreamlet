


#' Beeswarm plot of effect sizes for each assay
#' 
#' Beeswarm plot of effect sizes for each assay, colored by sign and FDR
#' 
#' @param res.dl \code{dreamletResult} object from \code{dreamlet()}
#' @param coef coefficient name fed to \code{topTable()}
#' @param fdr.range range for coloring FDR
#' @param assays which assays to plot
#' 
#' @return \code{ggplot2} of logFC by assay
#'
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce,
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#'
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~group_id)
#'
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data
#' res.dl <- dreamlet(res.proc, ~group_id)
#'
#' # Beeswarm plot of effect sizes for each assay, 
#' # colored by sign and FDR
#' plotBeeswarm( res.dl, 'group_idstim')
#
#' @export
#' @import ggplot2 
#' @importFrom dplyr arrange 
#' @importFrom ggbeeswarm geom_quasirandom
plotBeeswarm = function(res.dl, coef, fdr.range = 4, assays = assayNames(res.dl)){

	stopifnot(is(res.dl, 'dreamletResult'))

	adj.P.Val = logFC = score = assay = NULL

	# get results
	tab = topTable(res.dl, coef=coef, number=Inf) %>% 
			as_tibble %>%
			mutate(score = -log10(adj.P.Val) * sign(logFC)) %>%
			mutate(score = pmax(-fdr.range, pmin(fdr.range, score))) %>%
			filter(assay %in% assays) %>%
			mutate(assay = factor(assay, levels=assays))

	# make plot
	tab %>% 
		arrange(abs(score)) %>%  
		ggplot(aes(assay, logFC, color=score)) + 
		geom_quasirandom(method="maxout", size=1) +
		scale_color_gradient2(low = "blue", mid="grey", high="red", name = bquote(-log[10]~FDR ~ sign(logFC))) +
		theme_classic() +
		geom_hline(yintercept=0, color="black", linetype="dashed") +
		theme(legend.position = "bottom", aspect.ratio=1) +
		coord_flip()
}
