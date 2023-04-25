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
#' @param alpha transparency of points
#' @param assays which assays to plot
#' @param ... other arguments
#'
#' @return Plot of mean-variance trend
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
#' # Show mean-variance trend from voom
#' plotVoom(res.proc)
#' 
#' # plot for first two cell types
#' plotVoom(res.proc[1:2])
#' 
#' @export
#' @docType methods
#' @rdname plotVoom-methods
setGeneric("plotVoom", 
  function(x, ncol=3, alpha=.5,...){

  standardGeneric("plotVoom")
})


#' @rdname plotVoom-methods
#' @aliases plotVoom,dreamletProcessedData,dreamletProcessedData-method
setMethod("plotVoom", "dreamletProcessedData",
  function(x, ncol=3, alpha=.5, assays = names(x)){

	# Pass R CMD check
	y = NULL

	# intersect preserving order from assays
	assays = intersect(assays, names(x))
	if( length(assays) == 0) stop("No valid assays selected")

	# get common range across all plots
	###################################
	df_range = lapply( assays, function(id){

		if( !is.null(x[[id]]$voom.xy) ){
			res = with(x[[id]]$voom.xy, data.frame(range(x),range(y), id=id))
		}else{
			res = NULL
		}
		res
	})
	df_range = do.call(rbind, df_range)

	if( is.null(df_range) ){
		stop("Voom was not run on this object")
	}

	xlim = range(df_range$range.x.)
	ylim = c(0, max(df_range$range.y.))

	xlab = bquote(log[2](counts + 0.5))
	ylab = bquote(sqrt(standard~deviation))

	# only included assays were voom succeeded
	validAssays = droplevels(factor(unique(df_range)$id, assays))

	# make data.frame of points
	df.list = lapply( validAssays, function(id){
		with(x[[id]]$voom.xy, data.frame(id, x,y))
	})
	df_points = do.call(rbind, df.list)
	df_points$id = droplevels(factor(df_points$id, assays))
	df_points = df_points[order(df_points$id),]

	# make data.frame of curves
	df.list = lapply( validAssays, function(id){
		with(x[[id]]$voom.line, data.frame(id, x,y))
	})	
	df_curve = do.call(rbind, df.list)

	ggplot(df_points, aes(x,y)) + 
		geom_point(size=0.1, alpha=alpha) + 
		theme_classic() + 
		theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5))+ 
		facet_wrap(~id, ncol=ncol) + 
		xlab(xlab) + 
		ylab(ylab) + 
		xlim(xlim) + 
		ylim(ylim) + 
		geom_line(data = df_curve, aes(x,y), color="red")
})



#' @rdname plotVoom-methods
#' @aliases plotVoom,list,list-method
setMethod("plotVoom", "EList",
  function(x, ncol=3, alpha=.5){

	# Pass R CMD check
	y = NULL

	# get common range across all plots
	###################################

	df_range = with(x$voom.xy, data.frame(range(x),range(y)))

	xlim = range(df_range$range.x.)
	ylim = c(0, max(df_range$range.y.))

	xlab = bquote(log[2](counts + 0.5))
	ylab = bquote(sqrt(standard~deviation))

	# make data.frame of points
	df_points = with(x$voom.xy, data.frame(x,y))

	# make data.frame of curves
	df_curve = with(x$voom.line, data.frame(x,y))

	ggplot(df_points, aes(x,y)) + 
		geom_point(size=0.1, alpha=alpha) + 
		theme_classic() + 
		theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + xlab(xlab) + 
		ylab(ylab) + 
		xlim(xlim) + 
		ylim(ylim) + 
		geom_line(data = df_curve, aes(x,y), color="red")
})





