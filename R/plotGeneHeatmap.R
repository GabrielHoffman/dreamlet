#' Heatmap of genes and assays
#'
# Heatmap of genes and assays
#'
#' @param x A \code{dreamletResult} object
#' @param coef column number or column name specifying which coefficient or contrast of the linear model is of interest.
#' @param coef array of genes to include in plot
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(x)}
#' @param zmax maximum z.std value
#' @param transpose (default: FALSE) Use `coord_flip()` to flip axies
#'
#' @return Heatmap plot for specified genes and assays
#'
#' @export
#' @docType methods
#' @rdname plotGeneHeatmap-methods
setGeneric(
  "plotGeneHeatmap",
  function(x, coef, genes, assays = assayNames(x), zmax = NULL, transpose = FALSE) {
    standardGeneric("plotGeneHeatmap")
  }
)




#' Heatmap of genes and assays
#'
#' Heatmap of genes and assays
#'
#' @param x A \code{dreamletResult} object
#' @param coef column number or column name specifying which coefficient or contrast of the linear model is of interest.
#' @param genes array of genes to include in plot
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(x)}
#' @param zmax maximum z.std value
#' @param transpose (default: FALSE) Use `coord_flip()` to flip axies
#'
#' @return Heatmap plot for specified genes and assays
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
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~group_id)
#'
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data
#' res.dl <- dreamlet(res.proc, ~group_id)
#'
#' # Heatmap for specified subset of genes
#' plotGeneHeatmap(res.dl, coef = "group_idstim", genes = rownames(pb)[1:15])
#'
#' @rdname plotGeneHeatmap-methods
#' @aliases plotGeneHeatmap,dreamletResult,dreamletResult-method
#' @importFrom tidyr complete
setMethod(
  "plotGeneHeatmap", "dreamletResult",
  function(x, coef, genes, assays = assayNames(x), zmax = NULL, transpose = FALSE) {
    # extract gene-level results
    tab <- topTable(x, coef = coef, number = Inf)
    tab <- tab[tab$ID %in% genes, c("assay", "ID", "z.std", "adj.P.Val"), drop = FALSE]
    tab$ID <- factor(tab$ID, levels = genes)

    if (nrow(tab) == 0) stop("No genes retained")

    tab <- as.data.frame(tab[tab$assay %in% assays, , drop = FALSE])
    tab <- droplevels(tab)
    tab$assay <- factor(tab$assay, assays)

    if (nrow(tab) == 0) stop("No assays retained")

    # pass R CMD check
    assay <- ID <- z.std <- adj.P.Val <- NULL

    # fill empty values with NA
    tab <- as.data.frame(complete(tab, assay, ID))

    if (is.null(zmax)) {
      zmax <- max(abs(tab$z.std), na.rm = TRUE)
    } else {
      if (max(abs(tab$z.std), na.rm = TRUE) > zmax) {
        tab$z.std = pmax(tab$z.std, -zmax)
        tab$z.std = pmin(tab$z.std, zmax)
        warning("z.std > zmax", immediate. = TRUE)
      }
    }

    ncol <- length(unique(tab$assay))
    nrow <- length(unique(tab$ID))
    aspect.ratio <- nrow / ncol

    fig <- ggplot(tab, aes(assay, ID, fill = z.std, label = ifelse(adj.P.Val < 0.05, "*", ""))) +
      geom_tile() +
      theme_classic() +
      geom_text(vjust = 1, hjust = 0.5) +
      scale_fill_gradient2("z-statistic", low = "blue", mid = "white", high = "red", limits = c(-zmax, zmax), na.value = "grey70") +
      ylab("Gene") +
      xlab("Cell type")

    if (transpose) {
      fig <- fig + theme(aspect.ratio = 1 / aspect.ratio, axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5)) +
        coord_flip()
    } else {
      fig <- fig + theme(aspect.ratio = aspect.ratio, axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))
    }

    fig
  }
)
