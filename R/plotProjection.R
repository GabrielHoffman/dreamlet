# Gabriel Hoffman
# Jan 25, 2023


#' Plot 2D projection 
#' 
#' Plot 2D projection (i.e. UMAP, tSNE) for millions of cells efficiently
#' 
#' @param sce \code{SingleCellExperiment}
#' @param type field in \code{reducedDims(sce)} to plot
#' @param annotation column in \code{colData(sce)} to annotate each cell
#' @param pointsize Radius of rasterized point. Use \code{0} for single pixels(fastest).
#' @param pixels Vector with X and Y resolution of the raster, default \code{c(512,512)}
#' @param legend.position legend.position: the position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param text show \code{annotation} as text. Default \code{TRUE}
#' @param order specify order of levels for \code{annotation}
#' 
#' @details Uses \code{scattermore::geom_scattermore()} to plot millions of points efficiently
#' 
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#' 
#' data(example_sce)
#' 
#' plotProjection(example_sce, "TSNE", 'cluster_id', 1)
#' @import dplyr ggplot2 SingleCellExperiment
#' @importFrom scattermore geom_scattermore
#' @export
plotProjection = function(sce, type, annotation, pointsize = 0, pixels = c(512, 512), legend.position='none', text=TRUE, order){

	stopifnot(is(sce, "SingleCellExperiment"))
	stopifnot(type %in% names(reducedDims(sce)))
	stopifnot(annotation %in% colnames(colData(sce)))

	# PASS R CMD check
	CellType = X1 = X2 = x = y = NULL

	# Combine coordinates and annotations
	df = data.frame(reducedDim(sce, type), CellType = sce[[annotation]])
	colnames(df)[1:2] = c("X1", "X2")

	if( ! missing(order) ){
		df$CellType = factor(df$CellType, order)
	}

	# location of text
	df_txt = df %>%
	  group_by(CellType) %>%
	  summarize(x = median(X1), y = median(X2))

	fig = ggplot(df, aes(X1, X2, color=CellType)) +
	  geom_scattermore(pointsize = pointsize, pixels=pixels) + 
	  theme_classic() +
	  theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), axis.text = element_blank(), axis.ticks = element_blank(), legend.position=legend.position) +
	  xlab("Dim 1") +
	  ylab("Dim 2")

	if( text ){
		fig = fig + geom_text(data=df_txt, aes(x,y, label=CellType), color="black") 
	}  
	fig
}





