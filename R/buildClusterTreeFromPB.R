# Gabriel Hoffman
# Nov 30, 2022

#' Hierarchical clustering on cell types from pseudobulk
#'
#' Perform hierarchical clustering on cell types from pseudobulk by aggregating read counts from each cell type.
#'
#' @param pb \code{SingleCellObject} storing pseudobulk for each cell type in in \code{assay()} field
#' @param method clustering method for \code{hclust()}
#' @param dist.method distance metric
#' @param assays which assays to include
#'
#' @return hierarchical clustering object of class \code{hclust}
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
#' # Hierarchical clustering of cell types
#' hcl <- buildClusterTreeFromPB(pb)
#'
#' plot(hcl)
#'
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom edgeR filterByExpr DGEList calcNormFactors cpm
#' @importFrom stats dist hclust
#' @importFrom MatrixGenerics rowSums2
#' @export
buildClusterTreeFromPB <- function(pb, method = c("complete", "ward.D", "single", "average", "mcquitty", "median", "centroid", "ward.D2"), dist.method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), assays = assayNames(pb)) {
  method <- match.arg(method)
  dist.method <- match.arg(dist.method)

  # Combine counts for each cell cluster into a single point
  geneCounts <- lapply(assays, function(CT) {
    rowSums2(assay(pb, CT))
  })
  names(geneCounts) <- assays
  geneCounts <- do.call(cbind, geneCounts)

  # compute log2 CPM for genes with sufficient expression
  keep <- filterByExpr(geneCounts, group = rep(1, ncol(geneCounts)))
  dge <- DGEList(geneCounts[keep, , drop = FALSE])
  dge <- calcNormFactors(dge)
  geneExpr <- edgeR::cpm(dge, log = TRUE)

  # evaluate distance between pairs of cell clusters
  d <- dist(t(geneExpr), method = dist.method)

  hclust(d, method = method)
}
