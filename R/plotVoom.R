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
#' @import ggplot2
#' @export
#' @docType methods
#' @rdname plotVoom-methods
setGeneric("plotVoom", 
  function(x, ncol=3){

  standardGeneric("plotVoom")
})


#' @rdname plotVoom-methods
#' @aliases plotVoom,dreamletProcessedData,dreamletProcessedData-method
setMethod("plotVoom", "dreamletProcessedData",
  function(x, ncol=3){

	# Pass R CMD check
	y = NULL

	# get common range across all plots
	###################################
	df_range = lapply( names(x), function(id){

		if( !is.null(x[[id]]$geneExpr$voom.xy) ){
			res = with(x[[id]]$geneExpr$voom.xy, data.frame(range(x),range(y), id=id))
		}else{
			res = NULL
		}
		res
	})
	df_range = do.call(rbind, df_range)

	xlim = range(df_range$range.x.)
	ylim = c(0, max(df_range$range.y.))

	xlab = bquote(log[2](counts + 0.5))
	ylab = bquote(sqrt(standard~deviation))

	# only included assays were voom succeeded
	validAssays = unique(df_range)$id

	# make data.frame of points
	df.list = lapply( validAssays, function(id){
		with(x[[id]]$geneExpr$voom.xy, data.frame(id, x,y))
	})
	df_points = do.call(rbind, df.list)

	# make data.frame of curves
	df.list = lapply( validAssays, function(id){
		with(x[[id]]$geneExpr$voom.line, data.frame(id, x,y))
	})
	df_curve = do.call(rbind, df.list)

	ggplot(df_points, aes(x,y)) + geom_point(size=0.1) + theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + facet_wrap(~id, ncol=ncol) + xlab(xlab) + ylab(ylab) + xlim(xlim) + ylim(ylim) + geom_line(data = df_curve, aes(x,y), color="red")
})



#' @rdname plotVoom-methods
#' @aliases plotVoom,list,list-method
setMethod("plotVoom", "list",
  function(x, ncol=3){

	# Pass R CMD check
	y = NULL

	# get common range across all plots
	###################################

	df_range = with(x$geneExpr$voom.xy, data.frame(range(x),range(y)))

	xlim = range(df_range$range.x.)
	ylim = c(0, max(df_range$range.y.))

	xlab = bquote(log[2](counts + 0.5))
	ylab = bquote(sqrt(standard~deviation))

	# make data.frame of points
	df_points = with(x$geneExpr$voom.xy, data.frame(x,y))

	# make data.frame of curves
	df_curve = with(x$geneExpr$voom.line, data.frame(x,y))

	ggplot(df_points, aes(x,y)) + geom_point(size=0.1) + theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + xlab(xlab) + ylab(ylab) + xlim(xlim) + ylim(ylim) + geom_line(data = df_curve, aes(x,y), color="red")
})





