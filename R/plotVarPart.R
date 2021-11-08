
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
#' @param ... additional arguments
#' 
#' @return Violin plot showing variance fractions
#' 
#' @examples
#'  
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
#' genes = vp$gene[2:4]
#' plotVarPart(vp)
#'
#' @rdname plotVarPart-methods
#' @aliases plotVarPart,DataFrame,DataFrame-method
#' @importFrom reshape2 melt
#' @importFrom S4Vectors as.data.frame
#' @import ggplot2
setMethod("plotVarPart", "DataFrame",
	function( obj, col=c(ggColorHue(base::ncol(obj)-3), "grey85"), label.angle=20, main="", ylab = "",  convertToPercent = TRUE,...){

	if( any(!c("assay", "gene") %in% colnames(obj)) ){
		stop("obj must have columns with names: 'assay', and 'gene'")
	}	
			
	df = melt( as.data.frame(obj), id.vars = c("assay", "gene"))

	if( label.angle == '') label.angle = 20

	# pass R CMD check
	variable = value = NULL

	ggplot(df, aes(variable, value)) + 
		geom_violin( scale="width", aes(fill = factor(variable))) + 
		ylab("Variance explained (%)") + xlab('') + ylim(0, 1) + theme_bw() + 
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
		facet_wrap(~assay, ncol=3) + 
		ggtitle(main)
})
