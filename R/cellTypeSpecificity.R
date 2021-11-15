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
#' 	assay = "counts",    
#' 	cluster_id = 'cluster_id', 
#' 	sample_id = 'sample_id',
#' 	verbose=FALSE)
#' 
#' # Compute cell type specificity of each gene
#' df = cellTypeSpecificity( pb)
#' 
#' # Violin plot of specificity score for each cell type
#' plotViolin(df)
#' 
#' # Barplot of 5 genes
#' plotPercentBars( df, genes = rownames(df)[1:5])
#' 
#' # Compute the maximum specificity score for each gene
#' scoreMax = apply(df, 1, max)
#' head(scoreMax)
#' 
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom S4Vectors DataFrame
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

	df = DataFrame(geneExpr.fract)
	colnames(df) = colnames(geneExpr.fract)

	new("cellSpecificityValues", df)
}


#' Class cellSpecificityValues 
#'
#' Class \code{cellSpecificityValues} cell type specificity valupes for each gene and cell type
#'
#' @name cellSpecificityValues-class
#' @rdname cellSpecificityValues-class
#' @exportClass cellSpecificityValues
setClass("cellSpecificityValues", contains="DataFrame")
















