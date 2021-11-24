

#' Plot heatmap
#' 
#' Plot heatmap
#' 
#' @param x fractions for each gene
#' @param genes name of genes to plot
#' @param color color of heatmap
#'
#' @return heatmap
#'  
#' @export
#' @docType methods
#' @rdname plotHeatmap-methods
setGeneric("plotHeatmap", 
  function(x, genes = rownames(x), color="darkblue"){

  standardGeneric("plotHeatmap")
})



#' @export
#' @importFrom reshape2 melt
#' @import ggplot2 
#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,cellSpecificityValues,cellSpecificityValues-method
setMethod("plotHeatmap", "cellSpecificityValues",
  function(x, genes = rownames(x), color="darkblue"){

	# subset based on specified genes
	x = x[rownames(x) %in% unique(genes),]	

	# pass R CMD check
	value = variable = NA

  df = data.frame(gene = rownames(x), x, check.names=FALSE)

	df_melt = reshape2::melt(df, id.vars="gene")

	df_melt$gene = factor(df_melt$gene, unique(genes))
	df_melt$variable = factor(df_melt$variable, colnames(x))

	ratio = nlevels(df_melt$gene) / nlevels(df_melt$variable)

	# heatmap of cell type specificity
	ggplot(df_melt, aes(variable, gene, fill=value)) + 
		geom_tile() +
		theme_classic() + 
		theme(aspect.ratio=ratio, 
		  plot.title = element_text(hjust = 0.5),
		  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
		scale_fill_gradient(name = "Fraction of\nexpression", low="white", high=color, limits=c(0, 1)) +
		xlab('') + ylab('') +
		ggtitle("Cell type specificity scores")
})



#' @export
#' @importFrom reshape2 melt
#' @import ggplot2 
#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,matrix,matrix-method
setMethod("plotHeatmap", "matrix",
  function(x, genes = rownames(x), color="darkblue"){

	# subset based on specified genes
	x = x[rownames(x) %in% unique(genes),]	

	# pass R CMD check
	value = variable = NA

  df = data.frame(gene = rownames(x), x, check.names=FALSE)

	df_melt = reshape2::melt(df, id.vars="gene")

	df_melt$gene = factor(df_melt$gene, unique(genes))
	df_melt$variable = factor(df_melt$variable, colnames(x))

	ratio = nlevels(df_melt$gene) / nlevels(df_melt$variable)

	# heatmap of cell type specificity
	ggplot(df_melt, aes(variable, gene, fill=value)) + 
		geom_tile() +
		theme_classic() + 
		theme(aspect.ratio=ratio, 
		  plot.title = element_text(hjust = 0.5),
		  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
		scale_fill_gradient(name = bquote(log[2]~CPM), low="white", high=color) +
		xlab('') + ylab('') 
})





