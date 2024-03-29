# Gabriel Hoffman
# Nov 8, 2021

#' Class cellSpecificityValues
#'
#' Class \code{cellSpecificityValues} cell type specificity values for each gene and cell type
#'
#' @name cellSpecificityValues-class
#' @rdname cellSpecificityValues-class
#' @exportClass cellSpecificityValues
#' @return none
setClass("cellSpecificityValues", contains = "DFrame")




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
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#'
#' # Compute cell type specificity of each gene
#' df <- cellTypeSpecificity(pb)
#'
#' # Violin plot of specificity scores for each cell type
#' # Dashed line indicates genes that are equally expressed
#' # across all cell types.  For K cell types, this is 1/K
#' plotViolin(df)
#'
#' # Compute the maximum specificity score for each gene
#' scoreMax <- apply(df, 1, max)
#' head(scoreMax)
#'
#' # For each cell type, get most specific gene
#' genes <- rownames(df)[apply(df, 2, which.max)]
#'
#' # Barplot of 5 genes
#' plotPercentBars(df, genes = genes)
#'
#' # heatmap of 5 genes that are most cell type specific
#' dreamlet::plotHeatmap(df, genes = genes)
#'
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom S4Vectors DataFrame
#' @importFrom MatrixGenerics rowSums2
#' @export
cellTypeSpecificity <- function(pb, ...) {
  if (!is(pb, "SingleCellExperiment")) {
    stop("pb must be of class 'SingleCellExperiment'")
  }

  # sum counts for each cell type
  geneExpr <- lapply(assayNames(pb), function(key) {
    # get expression for this assay, and sum across all samples
    rowSums2(assay(pb, key), useNames = TRUE)
  })
  names(geneExpr) <- assayNames(pb)
  geneExpr <- do.call(cbind, geneExpr)
  rownames(geneExpr) <- rownames(pb)

  # identify genes with no reads
  idx <- which(rowSums2(geneExpr, useNames = TRUE) == 0)

  if (length(idx) > 0) {
    geneExpr <- geneExpr[-idx, , drop = FALSE]
    txt <- paste("There are", length(idx), "genes with no reads in this dataset. They are excluded here")
    warning(txt)
  }

  # get total expression counts for each cell,
  # and peform normalization
  dge <- DGEList(geneExpr)
  dge <- calcNormFactors(dge, ...)

  # evaluate counts per million
  geneExpr <- edgeR::cpm(dge, log = FALSE)

  df <- DataFrame(
    totalCPM = rowSums2(geneExpr, useNames = TRUE),
    geneExpr / rowSums2(geneExpr, useNames = TRUE), check.names = FALSE
  )

  new("cellSpecificityValues", df)
}
