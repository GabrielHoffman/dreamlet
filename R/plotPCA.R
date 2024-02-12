

setGeneric("plotPCA", BiocGenerics::plotPCA)



#' Plot PCA of gene expression for an assay
#' 
#' Compute PCA of gene expression for an assay, and plot samples coloring by outlier score
#' 
#' @param object \code{dreamletProcessedData} from \code{processAssays()}
#' @param assays assays / cell types to analyze
#' @param nPC number of PCs to uses for outlier score with \code{outlier()}
#' @param robust use robust covariance method, defaults to \code{FALSE}
#' @param ... arguments passed to \code{MASS::cov.rob()}
#' @param maxOutlierZ cap outlier z-scores at this value for plotting to maintain consistent color scale
#' @param nrow number of rows in plot
#' @param size size passed to \code{geom_point()}
#' 
#' @name plotPCA
#' @rdname plotPCA
#' @aliases plotPCA plotPCA,dreamletProcessedData-method
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
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~group_id)
#' 
#' # PCA to identify outliers
#' plotPCA( res.proc, "CD14+ Monocytes")
#
#' @seealso \code{outlierByAssay()}
#' @export
setMethod("plotPCA", signature(object="dreamletProcessedData"), function(object, assays = assayNames(object), nPC=2, robust = FALSE, ..., maxOutlierZ=20, nrow=2, size=2){

  stopifnot(all(assays %in% assayNames(object)))

  PC1 <- PC2 <- z <- pValue <- FDR <- NULL

  outlierByAssay( object, assays, robust = robust, ...) %>%
    tibble(FDR = p.adjust(pValue, "fdr")) %>%
    arrange(FDR) %>%
    ggplot(aes(PC1, PC2, 
      color = pmin(z, maxOutlierZ),
      shape = FDR < 0.05)) +
      geom_point(size=size) +
      theme_classic() +
      theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) +
      scale_color_gradient("Outlier z", limits=c(0, maxOutlierZ), low="grey60", high="red") +
      facet_wrap( ~ assay, nrow=nrow, scales="free")
})






