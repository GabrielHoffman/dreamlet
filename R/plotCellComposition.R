# Gabriel Hoffman
# Sept 28, 2021






#' Bar plot of cell compositions
#'
#' Bar plot of cell compositions
#'
#' @param obj matrix of [cells] x [samples] or \code{SingleCellExperiment} from \code{aggregateToPseudoBulk}
#' @param col array of colors.  If missing, use default colors.  If \code{names(col)} is the same as \code{arrayNames(obj)}, then colors will be assigned by assay name#'
#' @param width specify width of bars
#'
#' @return Barplot showing cell fractions
#'
#' @examples
#'
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
#' # show cell composition bar plots
#' plotCellComposition(pb)
#'
#' # extract cell counts
#' df_cellCounts <- cellCounts(pb)
#'
#' # show cell composition bar plots
#' plotCellComposition(df_cellCounts)
#'
#' @rdname plotCellComposition
#' @importFrom variancePartition plotPercentBars ggColorHue
#' @export
setGeneric(
  "plotCellComposition",
  function(obj, col, width = NULL) {
    standardGeneric("plotCellComposition")
  }
)



#' @export
#' @rdname plotCellComposition
#' @aliases plotCellComposition,SingleCellExperiment-method
setMethod(
  "plotCellComposition", "SingleCellExperiment",
  function(obj, col, width = NULL) {
    # extract cell counts and other meta-data
    df_cellCount <- do.call(rbind, int_colData(obj)$n_cells)

    .plotCellComposition(df_cellCount, col, width)
  }
)


#' @export
#' @rdname plotCellComposition
#' @aliases plotCellComposition,matrix-method
setMethod(
  "plotCellComposition", "matrix",
  function(obj, col, width = NULL) {
    .plotCellComposition(obj, col, width)
  }
)


#' @export
#' @rdname plotCellComposition
#' @aliases plotCellComposition,data.frame-method
setMethod(
  "plotCellComposition", "data.frame",
  function(obj, col, width = NULL) {
    .plotCellComposition(obj, col, width)
  }
)


#' @importFrom variancePartition plotPercentBars ggColorHue
.plotCellComposition <- function(countMatrix, col, width = NULL) {
  # compute fractions from counts
  df_frac <- apply(countMatrix, 1, function(x) x / sum(x))

  df <- as.data.frame(t(df_frac))

  if (missing(col)) {
    col <- ggColorHue(ncol(df))
  } else if (identical(sort(names(col)), sort(colnames(df)))) {
    col <- col[colnames(df)]
  } else if (length(col) < ncol(df)) {
    stop("Too few colors specified: ", length(col), " < ", ncol(df))
  }

  plotPercentBars(df, col = col, width = width) + ylab("Cell percentage")
}
