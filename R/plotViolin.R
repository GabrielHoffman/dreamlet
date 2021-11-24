

#' Plot Violins
#' 
#' Plot Violins
#' 
#' @param x fractions for each gene
#'
#' @return Violin plot
#'  
#' @export
#' @docType methods
#' @rdname plotViolin-methods
setGeneric("plotViolin", 
  function(x){

  standardGeneric("plotViolin")
})



#' @importFrom reshape2 melt
#' @import ggplot2 
#' @rdname plotViolin-methods
#' @aliases plotViolin,cellSpecificityValues,cellSpecificityValues-method
setMethod("plotViolin", "cellSpecificityValues",
  function(x){

	# pass R CMD check
	gene = value = variable = NA

	df = data.frame(gene = rownames(x), x, check.names=FALSE)

	df_melt = reshape2::melt(df, id.vars="gene")

	# Violin plot of cell type specificity
	ggplot(df_melt, aes(variable, value, fill=variable)) + 
		geom_violin() + 
		geom_boxplot(width=.1, fill="grey") + 
		theme_classic() + 
		theme(aspect.ratio=1, legend.position="none", 
		  plot.title = element_text(hjust = 0.5),
		  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
		scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
		ylab("Fraction of gene expression") + 
		 geom_hline(yintercept=1/ncol(df), color="grey60", linetype="dashed") +
		ggtitle("Cell type specificity scores") +
		xlab('')
})




