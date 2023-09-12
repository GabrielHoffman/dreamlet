#' Plot heatmap
#'
#' Plot heatmap
#'
#' @param x fractions for each gene
#' @param genes name of genes to plot
#' @param color color of heatmap
#' @param assays array of assays to plot
#' @param ... other arguments
#'
#' @return heatmap
#'
#' @export
#' @docType methods
#' @rdname plotHeatmap-methods
setGeneric(
  "plotHeatmap",
  function(x, genes = rownames(x), color = "darkblue", ...) {
    standardGeneric("plotHeatmap")
  }
)


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
#' # For each cell type, get most specific gene
#' genes <- rownames(df)[apply(df, 2, which.max)]
#'
#' # heatmap of 5 genes that are most cell type specific
#' dreamlet::plotHeatmap(df, genes = genes)
#' @export
#' @importFrom reshape2 melt
#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,cellSpecificityValues,cellSpecificityValues-method
setMethod(
  "plotHeatmap", "cellSpecificityValues",
  function(x, genes = rownames(x), color = "darkblue", assays = colnames(x)) {
    fig <- dreamlet::plotHeatmap(as.matrix(x)[, -1], genes, color, assays)

    fig +
      ggtitle("Cell type specificity scores") +
      scale_fill_gradient(name = "Fraction of\nexpression", low = "white", high = color, limits = c(0, 1))
  }
)


#' @export
#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,data.frame,data.frame-method
setMethod(
  "plotHeatmap", "data.frame",
  function(x, genes = rownames(x), color = "darkblue", assays = colnames(x)) {
    dreamlet::plotHeatmap(as.matrix(x), genes, color, assays)
  }
)



#' @export
#' @importFrom reshape2 melt
#' @rdname plotHeatmap-methods
#' @aliases plotHeatmap,matrix,matrix-method
setMethod(
  "plotHeatmap", "matrix",
  function(x, genes = rownames(x), color = "darkblue", assays = colnames(x)) {
    genes <- genes[!is.na(genes)]

    # intersect preserving order from assays
    assays <- intersect(assays, colnames(x))
    if (length(assays) == 0) stop("No valid assays selected")

    x <- x[, assays, drop = FALSE]

    # subset based on specified genes
    x <- x[rownames(x) %in% unique(genes), , drop = FALSE]

    # pass R CMD check
    value <- variable <- gene <- NA

    df <- data.frame(gene = rownames(x), x, check.names = FALSE)

    df_melt <- reshape2::melt(df, id.vars = "gene")

    df_melt$gene <- factor(df_melt$gene, unique(genes))
    df_melt$variable <- factor(df_melt$variable, assays)
    df_melt <- droplevels(df_melt)

    ratio <- nlevels(df_melt$gene) / nlevels(df_melt$variable)

    # heatmap of cell type specificity
    ggplot(df_melt, aes(variable, gene, fill = value)) +
      geom_tile() +
      theme_classic() +
      theme(
        aspect.ratio = ratio,
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)
      ) +
      scale_fill_gradient(name = "value", low = "white", high = color) +
      xlab("") +
      ylab("")
  }
)
