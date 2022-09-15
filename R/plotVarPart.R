
#' Violin plot of variance fractions
#'
#' Violin plot of variance fraction for each gene and each variable
#'
#' @param obj \code{varParFrac} object returned by \code{fitExtractVarPart} or \code{extractVarPart}
#' @param col vector of colors
#' @param label.angle angle of labels on x-axis
#' @param main title of plot
#' @param ylab text on y-axis
#' @param convertToPercent multiply fractions by 100 to convert to percent values
#' @param ncol number of columns in the plot
#' @param ... additional arguments
#' 
#' @return Violin plot showing variance fractions
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
#' # variance partitioning analysis
#' vp = fitVarPart( res.proc, ~ group_id)
#' 
#' # Summarize variance fractions genome-wide for each cell type
#' plotVarPart(vp)
#'
#' @rdname plotVarPart-methods
#' @aliases plotVarPart,DataFrame,DataFrame-method
#' @importFrom reshape2 melt
#' @importFrom S4Vectors as.data.frame
#' @import ggplot2
setMethod("plotVarPart", "DataFrame",
	function( obj, col=c(ggColorHue(base::ncol(obj)-3), "grey85"), label.angle=20, main="", ylab = "", convertToPercent = TRUE, ncol = 3,...){

	if( any(!c("assay", "gene") %in% colnames(obj)) ){
		stop("obj must have columns with names: 'assay', and 'gene'")
	}	
	
	# get assays when it is not defined in generic function		
	args <- list(...)
	if( 'assays' %in% names(args) ){
		assays = args$assays
	}else{
		assays = assayNames(obj)
	}
						
	df = melt( as.data.frame(obj[obj$assay %in% assays,]), id.vars = c("assay", "gene"))

	# if( label.angle == '') label.angle = 20

	# pass R CMD check
	variable = value = NULL

	fig = ggplot(df, aes(variable, 100*value)) + 
		geom_violin( scale="width", aes(fill = factor(variable)), color=NA) + 
		ylab("Variance explained (%)") + 
		xlab('') + 
		scale_y_continuous(limits=c(0, 100), expand=c(0,3)) +
		theme_classic() + 
		geom_boxplot(width=0.07, fill="grey", outlier.colour='black') + 
		scale_fill_manual(values=col) +
		theme(legend.position="none", plot.title=element_text(hjust=0.5),
			axis.text.x = element_text(size  = 13,
	                            angle = label.angle,
	                            hjust = 1,
	                            vjust = 1), 
			aspect.ratio=1,
			text 	= element_text(colour="black"), 
			axis.text 	= element_text(colour="black"),
			legend.text = element_text(colour="black")) + 
		ggtitle(main)

	if( length(unique(df$assay)) > 1 ){
		fig = fig + facet_wrap(~assay, ncol=ncol)
	}

	fig 
})
