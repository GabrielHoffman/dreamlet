

setGeneric("plotPCA", BiocGenerics::plotPCA)



#' Plot PCA of gene expression for an assay
#' 
#' Compute PCA of gene expression for an assay, and plot samples coloring by outlier score
#' 
#' @param object \code{dreamletProcessedData} from \code{processAssays()}
#' @param assays assays / cell types to analyze
#' @param nPC number of PCs to uses for outlier score with \code{outlier()}
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
#' @export
setMethod("plotPCA", signature(object="dreamletProcessedData"), function(object, assays = assayNames(object), nPC=2, maxOutlierZ=20, nrow=2, size=1){

  stopifnot(all(assays %in% assayNames(object)))

  PC1 <- PC2 <- z <- NULL

  df = lapply(assays, function(id){
    # get normalized expression
    Y <- assay(object, id)$E

    # PCA
    dcmp <- prcomp(scale(t(Y)))

    # outlier analysis on first 2 PCs
    df_outlier <- outlier(dcmp$x[,seq(nPC)] * dcmp$sdev[seq(nPC)], FALSE)

    # plot
    dcmp$x[,seq(2)] %>%
      data.frame(df_outlier, assay = id) %>%
      as_tibble 
  })
  df = bind_rows(df)

   df %>%
    arrange(z) %>%
    ggplot(aes(PC1, PC2, color = pmin(z, maxOutlierZ))) +
      geom_point(size=size) +
      theme_classic() +
      # ggtitle(assays) + 
      theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) +
      scale_color_gradient("Outlier z", limits=c(0, maxOutlierZ), low="grey40", high="red") +
      facet_wrap(~assay, nrow=nrow)
})


