#' Class vpDF
#'
#' Class \code{vpDF} stores results for each gene for each assay
#'
#' @name vpDF-class
#' @rdname vpDF-class
#' @exportClass vpDF
#' @importFrom S4Vectors DataFrame
#' @return none
setClass("vpDF", contains = "DFrame", slots = c(df_details = "data.frame", errors = "list", error.initial = "list"))



#' Get assayNames
#'
#' Get assayNames
#'
#' @param x vpDF object
#' @param ... additional arguments
#'
#' @rdname assayNames-methods
#' @aliases assayNames,vpDF,vpDF-method
#' @export
setMethod(
  "assayNames", signature(x = "vpDF"),
  function(x, ...) {
    levels(x$assay)
  }
)

#' Get assays by name
#'
#' Get assays by name
#'
#' @param x vpDF object
#' @param i number indicating index, or string indicating assay
#' @param withDimnames not used
#'
#' @rdname assay-methods
#' @aliases assay,vpDF,vpDF-method
#' @export
setMethod(
  "assay", signature(x = "vpDF"),
  function(x, i, withDimnames = TRUE, ...) {
    if (is.numeric(i)) {
      i <- assayNames(x)[i]
    }
    x[x$assay == i, ]
  }
)

#' Sort variance partition statistics
#'
#' Sort variance partition statistics
#'
#' @param x object returned by \code{fitVarPart()}
#' @param FUN function giving summary statistic to sort by.  Defaults to sum
#' @param decreasing  logical.  Should the sorting be increasing or decreasing?
#' @param last columns to be placed on the right, regardless of values in these columns
#' @param ... other arguments to sort
#'
#' @return \code{data.frame} with columns sorted by mean value, with Residuals in last column
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
#' # variance partitioning analysis
#' vp <- fitVarPart(res.proc, ~group_id)
#'
#' # Summarize variance fractions genome-wide for each cell type
#' plotVarPart(sortCols(vp))
#'
#' @importMethodsFrom variancePartition sortCols
#' @importFrom stats median
#' @export
#' @rdname sortCols-method
#' @aliases sortCols,vpDF-method
setMethod(
  "sortCols", "vpDF",
  function(x, FUN = sum, decreasing = TRUE, last = c("Residuals", "Measurement.error"), ...) {
    if (nrow(x) == 0) {
      stop("vpDF object has no rows")
    }

    # perform storting without the first two annotation columns
    res <- sortCols(as.data.frame(x[, -c(1, 2), drop = FALSE]), FUN, decreasing, last, ...)

    # add the annotation columns back to the sorted data.frame
    new("vpDF", DataFrame(x[, c(1, 2)], res), df_details = x@df_details)
  }
)



#' @export
#' @rdname details-methods
#' @aliases details,vpDF-method
setMethod(
  "details", "vpDF",
  function(object) {
    object@df_details
  }
)

#' @export
#' @rdname seeErrors-methods
#' @aliases seeErrors,vpDF-method
#' @importFrom dplyr as_tibble
setMethod(
  "seeErrors", "vpDF",
  function(obj, initial = FALSE) {
    if (!initial) {
      df <- lapply(names(obj@errors), function(id) {
        if (length(obj@errors[[id]]) == 0) {
          return(NULL)
        }
        data.frame(
          assay = id,
          feature = names(obj@errors[[id]]),
          errorText = obj@errors[[id]]
        )
      })
    } else {
      df <- lapply(names(obj@error.initial), function(id) {
        message(id)
        if (length(obj@error.initial[[id]]) == 0) {
          return(NULL)
        }
        data.frame(
          assay = id,
          errorTextInitial = obj@error.initial[[id]]
        )
      })
    }
    as_tibble(do.call(rbind, df))
  }
)
