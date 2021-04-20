# Gabriel Hoffman
# April 6, 2021
#
# plot voom applied to pseudo-bulk from each cell type


#' Plot voom curves from each cell type
#'
#' Plot voom curves from each cell type
#' 
#' @param x dreamletProcessedData
#' @param ncol number of columns in the plot
#'
#' @import ggplot2 cowplot
#' @export
plotVoom = function( x, ncol=3){

	stopifnot(is(x, "dreamletProcessedData"))

	# Pass R CMD check
	y = NULL

	# get common range across all plots
	###################################
	df_range = lapply( names(x), function(id){

		with(x[[id]]$geneExpr$voom.xy, data.frame(range(x),range(y)))
	})
	df_range = do.call(rbind, df_range)

	xlim = range(df_range$range.x.)
	ylim = range(df_range$range.y.)

	xlab = bquote(log[2](counts + 0.5))
	ylab = bquote(sqrt(standard~deviation))

	# make list of plots
	figList = lapply( names(x), function(id){

		df = with(x[[id]]$geneExpr$voom.xy, data.frame(x,y))
		df_curve = with(x[[id]]$geneExpr$voom.line, data.frame(x,y))

		ggplot(df, aes(x, y)) + geom_point(size=0.1) + geom_line(data=df_curve, aes(x,y), color="red") + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + xlab(xlab) + ylab(ylab) + ggtitle(id) + xlim(xlim) + ylim(ylim)
	})

	plot_grid(plotlist = figList, ncols=ncol, align="hv", axis='tblr')
}