# extract table of cell counts from 'int_colData'
# of pseudobulks as returned by 'aggregateData'
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment int_colData
.n_cells <- function(x) {
  y <- int_colData(x)$n_cells
  if (is.null(y)) {
    return(NULL)
  }
  if (length(metadata(x)$agg_pars$by) == 2) {
    y <- as.matrix(data.frame(y, check.names = FALSE))
  }
  return(as.table(y))
}


#' Check variables in a formula
#'
#' Check that variables in formula are present in the data
#'
#' @param formula formula of variables to check
#' @param data data.frame storing variables in the formula
#'
#' @return If formula is valid, return TRUE.  Else throw error
#'
#' @examples
#'
#' # Valid formula
#' dreamlet:::checkFormula(~speed, cars)
#'
#' # Not valid formula
#' # dreamlet:::checkFormula( ~ speed + a, cars)
#'
#' @importFrom stats terms
checkFormula <- function(formula, data) {
  stopifnot(is(formula, "formula"))
  stopifnot(is(data, "data.frame"))

  v <- all.vars(formula)
  found <- v %in% colnames(data)

  if (any(!found)) {
    txt <- paste("Variables in formula are not found in data:\n   ", paste(v[!found], collapse = ", "))
    stop(txt)
  }
}


# Check if formula is full rank
#' @importFrom lme4 nobars
isFullRank <- function(formula, data) {
  design <- model.matrix(nobars(formula), data)

  qr(design)$rank >= ncol(design)
}




#' Check if two formulas are equal
#'
#' Check if two formulas are equal by evaluating the formulas and extracting terms
#'
#' @param formula1 first formula
#' @param formula2 second formula
#'
#' @return boolean value indciating of formulas are equivalent
#'
#' @examples
#'
#' # These formulas are equivalent
#' formula1 <- ~ Size + 1
#' formula2 <- ~ 1 + Size
#'
#' dreamlet:::equalFormulas(formula1, formula2)
#'
#' @importFrom stats terms
equalFormulas <- function(formula1, formula2) {
  # extract terms from forumula1
  trmf1 <- terms(as.formula(formula1))
  intercept1 <- attr(trmf1, "intercept")
  fterms1 <- attr(trmf1, "term.labels")

  # extract terms from forumula2
  trmf2 <- terms(as.formula(formula2))
  intercept2 <- attr(trmf2, "intercept")
  fterms2 <- attr(trmf2, "term.labels")

  # check equality of intercept and variables
  # sort because the order might be different
  (intercept1 == intercept2) & identical(sort(fterms1), sort(fterms2))
}


#' Extract residuals from \code{dreamletResult}
#'
#' Extract residuals from \code{dreamletResult}
#'
#' @param object \code{dreamletResult} object
#' @param y \code{dreamletProcessedData} object
#' @param ... other arguments
#' @param type compute either \code{"response"} residuals or \code{"pearson"} residuals.
#'
#' @details \code{"response"} residuals are the typical residuals returned from \code{lm()}. \code{"pearson"} residuals divides each residual value by its estimated standard error.  This requires specifying \code{y}
#'
#' @return residuals from model fit
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
#' # extract typical residuals for each assay (i.e. cell type)
#' # Return list with entry for each assay with for retained samples and genes
#' resid.lst <- residuals(res.dl)
#'
#' # Get Pearson residuals:
#' # typical residuals scaled by the standard deviation
#' residPearson.lst <- residuals(res.dl, res.proc, type = "pearson")
#'
#' @importMethodsFrom BiocGenerics residuals
#' @rdname residuals-methods
#' @aliases residuals,dreamletResult,dreamletResult-method
#' @export
setMethod(
  "residuals", "dreamletResult",
  function(object, y, ..., type = c("response", "pearson")) {
    type <- match.arg(type)
    hasy <- !missing(y)

    if (hasy) {
      y.type <- is(y, "dreamletProcessedData")
    } else {
      y.type <- FALSE
    }

    if (type == "pearson" & (!hasy | !y.type)) {
      stop("When type is 'pearson' specify y as dreamletProcessedData object")
    }

    if (type == "pearson") {
      if (!all(assayNames(object) %in% assayNames(y))) {
        stop("Not all cell types object are present in y")
      }
    }

    res <- lapply(names(object), function(key) {
      if (hasy && type == "pearson") {
        residuals(object[[key]], y = y[[key]], ..., type = type)
      } else {
        residuals(object[[key]], ..., type = type)
      }
    })
    names(res) <- names(object)

    res
  }
)


#' @importFrom dplyr bind_rows
get_metadata_aggr_means <- function(x) {
  # when SingleCellExperiments are cbind'ed,
  # the entries in metadata(x) are added as a list instead of rbind'ed
  # here, rbind the aggr_means lists entries
  idx <- which(names(metadata(x)) == "aggr_means")

  bind_rows(metadata(x)[idx])
}
