# Gabriel Hoffman
# Nov 8, 2021



#' Get cell type specificity of gene expression
#'
#' For each gene, compute fraction of overall expression attributable to each cell type
#' 
#' @param pb \code{SingleCellExperiment} of pseudobulk data where easy \code{assay} is a cell type.
#' @param ... other arguments passed to \code{edgeR::calcNormFactors()}
#' 
#' @return matrix of the fraction of expression attributable to each cell type for each gene.
#'
#' @details Sum counts for each cell type, and compute the fraction of counts-per-million attributable to each cell type for each gene
#' 
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#' 
#' data(example_sce)
#' 
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce, 
#' assay = "counts",    
#' cluster_id = 'cluster_id', 
#' sample_id = 'sample_id',
#' verbose=FALSE)
#' 
#' # Compute cell type specificity of each gene
#' df = cellTypeSpecificity( pb)
#' 
#' df_melt = reshape2::melt(df)
#' 
#' # Violin plot of cell type specificity
#' ggplot(df_melt, aes(Var2, value, fill=Var2)) + 
#'	geom_violin() + 
#'	geom_boxplot(width=.1, fill="grey") + 
#'	theme_classic() + 
#'	theme(aspect.ratio=1, legend.position="none", 
#'		plot.title = element_text(hjust = 0.5),
#'		axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
#'	ylim(0,1) + 
#'	ylab("Fraction of gene expression") + 
#'	ggtitle("Cell type specificity")
#' 
#' # Barplot of 5 genes
#' genes = rownames(df)[1:5]
#' 
#' ggplot(df_melt[df_melt$Var1 %in% genes,], aes(Var1, value, fill=Var2)) + 
#' 	geom_bar(stat="identity") + 
#' 	theme_classic() + 
#' 	theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + 
#' 	coord_flip() + 
#' 	scale_fill_discrete(name = "Cell type") + 
#' 	xlab("Gene") + 
#' 	ylab("Fraction of gene expression") + 
#' 	scale_y_continuous(limits=c(0, 1), expand=c(0,0)) + 
#'  ggtitle("Cell type specificity")
#' 
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @export
cellTypeSpecificity = function(pb,...){

	# sum counts for each cell type
	geneExpr = lapply( assayNames(pb), function(key){

		# get expression for this assay, and sum across all samples
		rowSums(assay(pb, key))
		})
	names(geneExpr) = assayNames(pb)
	geneExpr = do.call(cbind, geneExpr)

	# get total expression counts for each cell,
	# and peform normalization
	dge = DGEList(geneExpr)
	dge = calcNormFactors(dge,...)

	# evalute counts per million
	geneExpr = cpm(dge, log=FALSE)

	# for each gene, compute fraction of expression
	geneExpr.fract = t(apply(geneExpr, 1, function(x) x/ sum(x)))

	geneExpr.fract
}



