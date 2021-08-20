



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
