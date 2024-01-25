#' Class dreamletProcessedData
#'
#' Class \code{dreamletProcessedData}
#'
#' @name dreamletProcessedData-class
#' @rdname dreamletProcessedData-class
#' @exportClass dreamletProcessedData
#' @return none
setClass("dreamletProcessedData", contains = "list", slots = c(data = "data.frame", metadata = "data.frame", by = "vector", df_details = "data.frame", errors = "list", error.initial = "list"))

#' Subset with brackets
#'
#' Subset with brackets
#'
#' @param x \code{dreamletProcessedData} object
#' @param i indeces to extract
#'
#' @return entries stored at specified index

#' @rdname extract-methods
#' @aliases [,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod(
  "[", signature(x = "dreamletProcessedData"),
  function(x, i) {
    res <- new("dreamletProcessedData", x@.Data[i],
      data = x@data,
      metadata = x@metadata,
      by = x@by
    )
    names(res) <- names(x)[i]
    res
  }
)


setGeneric("assayNames", SummarizedExperiment::assayNames)
setGeneric("assay", SummarizedExperiment::assay)
setGeneric("colData", SummarizedExperiment::colData)
setGeneric("colData<-", SummarizedExperiment::`colData<-`)
setGeneric("metadata", S4Vectors::metadata)

#' Get assayNames
#'
#' Get assayNames
#'
#' @param x \code{dreamletProcessedData} object
#' @param ... other arguments
#'
#' @return array of assay names
#' @rdname assayNames-methods
#' @aliases assayNames,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod(
  "assayNames", signature(x = "dreamletProcessedData"),
  function(x, ...) {
    names(x)
  }
)

#' Get assay
#'
#' Get assay
#'
#' @param x \code{dreamletProcessedData} object
#' @param i number indicating index, or string indicating assay
#' @param withDimnames not used
#' @param ... other arguments
#'
#' @return return ith assay
#'
#' @rdname assay-methods
#' @aliases assay,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod(
  "assay", signature(x = "dreamletProcessedData"),
  function(x, i, withDimnames = TRUE, ...) {
    x[[i]]
  }
)


#' Extract colData from \code{dreamletProcessedData}
#'
#' Extract colData from \code{dreamletProcessedData}
#'
#' @param x A \code{dreamletProcessedData} object
#' @param ... other arguments
#'
#' @return object from \code{colData} field
#' @export
setMethod(
  "colData", "dreamletProcessedData",
  function(x, ...) {
    x@data
  }
)



#' Set colData
#'
#' Set colData of dreamletProcessedData, and check for same dimensions and rownames
#'
#' @param x \code{dreamletProcessedData} object
#' @param ... other arguments
#' @param value \code{data.frame} or object that can be coerced to it
#'
#' @return none
#' @export
setMethod(
  "colData<-", "dreamletProcessedData",
  function(x, ..., value) {
    # convert to data.frame
    value <- as.data.frame(value)

    # check dimensions
    if (nrow(x@data) != nrow(value)) {
      stop("Number of rows in colData(x) must remain the same")
    }

    # check same rownames
    if (!all.equal(rownames(x@data), rownames(value))) {
      stop("rownames(colData(x)) must remain the same in the new data.frame")
    }

    x@data <- value
    x
  }
)

#' Extract metadata from \code{dreamletProcessedData}
#'
#' Extract metadata from \code{dreamletProcessedData}
#'
#' @param x A dreamletProcessedData object
#'
#' @return object from \code{metadata} field
#' @aliases metadata,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod(
  "metadata", "dreamletProcessedData",
  function(x) {
    x@metadata
  }
)


#' Show object
#'
#' Show object
#'
#' @param object \code{dreamletProcessedData} object
#'
#' @rdname show-methods
#' @importFrom utils head tail
#' @importFrom S4Vectors coolcat
#' @aliases show,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod(
  "show", "dreamletProcessedData",
  function(object) {
    cat("class:", class(object), "\n")

    # assay
    nms <- assayNames(object)
    if (is.null(nms)) {
      nms <- character(length(assays(object, withDimnames = FALSE)))
    }
    coolcat("assays(%d): %s\n", nms)

    # colData
    nms <- names(colData(object))
    if (is.null(nms)) {
      nms <- character(length(colData(object, withDimnames = FALSE)))
    }
    coolcat("colData(%d): %s\n", nms)

    # metadata
    nms <- names(metadata(object))
    if (is.null(nms)) {
      nms <- character(length(metadata(object, withDimnames = FALSE)))
    }
    coolcat("metadata(%d): %s\n", nms)

    df_count <- lapply(object, function(obj) dim(obj))
    df_count <- do.call(rbind, df_count)

    if (is.null(df_count)) {
      cat("No assays retained\n")
    } else {
      cat("Samples:\n min:", min(df_count[, 2]), "\n max:", max(df_count[, 2]))
      cat("\nGenes:\n min:", min(df_count[, 1]), "\n max:", max(df_count[, 1]), "\n")

      # metadata
      nms <- names(details(object))
      if (is.null(nms)) {
        nms <- character(length(metadata(object, withDimnames = FALSE)))
      }
      coolcat("details(%d): %s\n", nms)
    }
  }
)



#' Print object
#'
#' Print object
#'
#' @param x \code{dreamletProcessedData} object
#' @param ... other arguments
#'
#' @export
#' @rdname print-methods
#' @aliases print,dreamletProcessedData,dreamletProcessedData-method
setMethod(
  "print", "dreamletProcessedData",
  function(x, ...) {
    show(x)
  }
)

# # setGene
# #' Extract a subset of samples
# #'
# #' Extract a subset of samples
# #'
# #' @param x dreamletProcessedData
# #' @param ids column names to retain
# #'
# #' @export
# subsetSamples = function(x, ids){

# 	stopifnot( is(x, 'dreamletProcessedData'))

# 	# for each assay
# 	for(i in seq_len(length(x)) ){

# 		# intersect ids with column names
# 		include = intersect(ids, colnames(x[[i]]))

# 		# extract samples with these column names
# 		x[[i]] = x[[i]][,include]
# 	}

# 	x
# }



#' Extract details from dreamletProcessedData
#'
#' Extract details from \code{dreamletProcessedData}
#'
#' @param object A \code{dreamletProcessedData} object
#'
#' @return Extract detailed information from some classes
#'
#' @importMethodsFrom GSEABase details
#' @rdname details-methods
#' @export
setGeneric("details", getGeneric("details", package = "GSEABase"))


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
#' # For each cell type, number of samples retained,
#' # and variables retained
#' details(res.proc)
#'
#' @export
#' @rdname details-methods
#' @aliases details,dreamletProcessedData-method
setMethod(
  "details", "dreamletProcessedData",
  function(object) {
    object@df_details
  }
)


#' @export
#' @rdname details-methods
#' @aliases details,dreamletResult-method
setMethod(
  "details", "dreamletResult",
  function(object) {
    object@df_details
  }
)





#' Extract normalized expression and \code{colData}
#'
#' Extract normalized expression and \code{colData}
#'
#' @param x \code{dreamletProcessedData} object
#' @param assay assay to extract
#' @param cols columns in \code{colData(x)} to extract.  defaults to all columns as \code{colnames(colData(x))}
#' @param genes genes to extract from \code{assay(x, assay)$E}. defaults to all genes as \code{rownames(x)}
#'
#' @rdname extractData-methods
#' @export
setGeneric("extractData", function(x, assay, cols = colnames(colData(x)), genes = rownames(x)) standardGeneric("extractData"))


#' Extract normalized expression and \code{colData}
#'
#' Extract normalized (i.e. log2 CPM) expression and \code{colData} from \code{dreamletProcessedData}
#'
#' @param x \code{dreamletProcessedData} object
#' @param assay assay to extract
#' @param cols columns in \code{colData(x)} to extract.  defaults to all columns as \code{colnames(colData(x))}
#' @param genes genes to extract from \code{assay(x, assay)$E}. defaults to all genes as \code{rownames(x)}
#'
#' @return \code{data.frame} or \code{DataFrame} of merged expression and colData
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
#' # Extract all:
#' # Extract tibble of colData merged with expression.
#' # variables and genes are stored as columns, samples as rows
#' df_merge <- extractData(res.proc, "B cells")
#'
#' # first few columns
#' df_merge[, 1:6]
#'
#' # Extract subset:
#' df_merge <- extractData(res.proc, "B cells", cols = "group_id", genes = c("SSU72", "U2AF1"))
#'
#' df_merge
#'
#' # Boxplot of expression
#' boxplot(SSU72 ~ group_id, df_merge)
#' #
#' @importFrom S4Vectors merge
#' @importFrom dplyr as_tibble
#' @rdname extractData-methods
#' @aliases extractData,dreamletProcessedData-method
#' @export
setMethod(
  "extractData", c(x = "dreamletProcessedData", assay = "character"),
  function(x, assay, cols = colnames(colData(x)), genes = rownames(assay(x, assay))) {
    # Check requested assay
    if (!assay %in% assayNames(x)) {
      stop("assay not found: ", assay)
    }

    # Check requested columns of colData(x)
    cols <- unique(cols)
    notFound <- cols[!cols %in% colnames(colData(x))]
    if (length(notFound) > 0) {
      txt <- paste(notFound[seq(min(5, length(notFound)))], collapse = ", ")
      stop("columns not found: ", txt)
    }

    # Check requested genes
    genes <- unique(genes)
    notFound <- genes[!genes %in% rownames(assay(x, assay))]
    if (length(notFound) > 0) {
      txt <- paste(notFound[seq(min(5, length(notFound)))], collapse = ", ")
      stop("genes not found: ", txt)
    }

    # merge data
    df <- merge(colData(x)[, cols, drop = FALSE],
      t(assay(x, assay)$E[genes, , drop = FALSE]),
      by = "row.names"
    )

    as_tibble(df)
  }
)

#' @export
#' @rdname seeErrors-methods
#' @aliases seeErrors,dreamletProcessedData-method
#' @importFrom dplyr as_tibble
setMethod(
  "seeErrors", "dreamletProcessedData",
  function(obj) {

    # Initial fit
    df <- lapply(names(obj@error.initial), function(id) {
      if (length(obj@error.initial[[id]]) == 0) {
        return(NULL)
      }
      tibble(
        assay = id,
        errorTextInitial = obj@error.initial[[id]]
      )
    })
    df <- bind_rows(df)

    txt = paste("   Assay-level errors:", nrow(df))
    message(txt)

    # Gene-level
    df2 <- lapply(names(obj@errors), function(id) {
      if (length(obj@errors[[id]]) == 0) {
        return(NULL)
      }
      tibble(
        assay = id,
        feature = names(obj@errors[[id]]),
        errorText = obj@errors[[id]]
      )
    })
    df2 <- bind_rows(df2)

    txt = paste("   Gene-level errors:", nrow(df2))
    message(txt)

    list(assayLevel = df, geneLevel = df2)
  }
)
