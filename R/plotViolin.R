

#' Plot Violins
#' 
#' Plot Violins
#' 
#' @param x fractions for each gene
#' @param assays array of assays to plot
#' @param ... other arguments
#'
#' @return Violin plot
#'  
#' @export
#' @docType methods
#' @rdname plotViolin-methods
setGeneric("plotViolin", 
  function(x,...){

  standardGeneric("plotViolin")
})


#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#' 
#' data(example_sce)
#' 
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce, 
#' 	assay = "counts",    
#' 	cluster_id = 'cluster_id', 
#' 	sample_id = 'sample_id',
#' 	verbose=FALSE)
#' 
#' # Compute cell type specificity of each gene
#' df = cellTypeSpecificity( pb)
#' 
#' # Violin plot of specificity scores for each cell type
#' # Dashed line indicates genes that are equally expressed
#' # across all cell types.  For K cell types, this is 1/K
#' plotViolin(df)
#' @importFrom reshape2 melt
#' @import ggplot2 
#' @rdname plotViolin-methods
#' @aliases plotViolin,cellSpecificityValues,cellSpecificityValues-method
setMethod("plotViolin", "cellSpecificityValues",
  function(x, assays=colnames(x)){

	# pass R CMD check
	gene = value = variable = NA

	# omit column totalCPM, if it exists
	i = which(colnames(x) == "totalCPM")
	if( length(i) > 0) x = x[,-1]

	# intersect preserving order from assays
	assays = intersect(assays, colnames(x))
	if( length(assays) == 0) stop("No valid assays selected")
	x = x[,assays,drop=FALSE]

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




