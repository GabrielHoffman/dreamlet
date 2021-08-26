



#' Class vpDF
#'
#' Class \code{vpDF} stores results for each gene for each assay
#'
#' @name vpDF-class
#' @rdname vpDF-class
#' @exportClass vpDF
setClass("vpDF", contains="DataFrame")

#' Get assayNames
#' 
#' Get aassayNames
#' 
#' @param x vpDF object
#' @param ... additional arguments
#'
#' @rdname vpDF-class
#' @export
setMethod("assayNames", signature(x="vpDF"),
	function(x, ...){   
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
#' @rdname vpDF-class
#' @export
setMethod("assay", signature(x="vpDF"),
	function(x, i, withDimnames=TRUE,...){ 
		if( is.numeric(i) ){
			i = assayNames(x)[i]
		}
		x[x$assay == i,]
	}
)

#' Sort variance partition statistics
#'
#'
#' @param x object returned by \code{extractVarPart()} or \code{fitExtractVarPartModel()}
#' @param FUN function giving summary statistic to sort by.  Defaults to median
#' @param decreasing  logical.  Should the sorting be increasing or decreasing?  
#' @param last columns to be placed on the right, regardless of values in these columns
#' @param ... other arguments to sort 
#'
#' @export
#' @rdname sortCols-method
#' @aliases sortCols,vpDF-method
#' @importFrom stats median
setMethod("sortCols", "vpDF",
	function( x, FUN=median, decreasing = TRUE, last=c("Residuals", "Measurement.error"), ... ){
 		
 		# perform storting without the first two annotation columns
		res = sortCols(as.data.frame(x[,-c(1,2)]), FUN, decreasing, last, ... )

		# add the annotation columns back to the sorted data.frame
		new("vpDF", DataFrame(x[,c(1,2)], res))
 	}
)


