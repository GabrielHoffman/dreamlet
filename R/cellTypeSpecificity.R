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
#' # heatmap of 5 genes that are most cell type specific
#' genes = rownames(df)[apply(df, 2, which.max)]
#' dreamlet::plotHeatmap( df, genes = genes)
#' 
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom S4Vectors DataFrame
#' @export
cellTypeSpecificity = function(pb,...){

	if( ! is(pb, 'SingleCellExperiment') ){
		stop("pb must be of class 'SingleCellExperiment'")
	}

	# sum counts for each cell type
	geneExpr = lapply( assayNames(pb), function(key){

		# get expression for this assay, and sum across all samples
		rowSums(assay(pb, key))
		})
	names(geneExpr) = assayNames(pb)
	geneExpr = do.call(cbind, geneExpr)

	# identify genes with no reads
	idx = which(rowSums(geneExpr) == 0)

	if(length(idx) > 0){

		geneExpr = geneExpr[-idx,]
		txt = paste("There are", length(idx), "genes with no reads in this dataset. They are excluded here")
		warning(txt)
	}

	# get total expression counts for each cell,
	# and peform normalization
	dge = DGEList(geneExpr)
	dge = calcNormFactors(dge,...)

	# evalute counts per million
	geneExpr = cpm(dge, log=FALSE)

	# for each gene, compute fraction of expression
	# geneExpr.fract = t(apply(geneExpr, 1, function(x) x/ sum(x)))
	# df = DataFrame(geneExpr.fract)
	# colnames(df) = colnames(geneExpr.fract)

	df = DataFrame(geneExpr / rowSums(geneExpr), check.names=FALSE)

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
















