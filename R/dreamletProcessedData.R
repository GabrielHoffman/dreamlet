
#' Class dreamletProcessedData
#'
#' Class \code{dreamletProcessedData} 
#'
#' @name dreamletProcessedData-class
#' @rdname dreamletProcessedData-class
#' @exportClass dreamletProcessedData
setClass("dreamletProcessedData", contains="list", slots = c(data = 'data.frame', metadata='data.frame', pkeys="vector"))

# Convert SingleCellExperiment to dreamletProcessedData
# @export
# setAs("SingleCellExperiment", "dreamletProcessedData", function(from){

# 	assayLists = lapply(assayNames(from), function(k){
# 		obj = new("EList", list(E=assay(from, k)))
# 		obj$isCounts = FALSE
# 		obj
# 	})
# 	names(assayLists) = assayNames(from)

# 	new("dreamletProcessedData", assayLists, data = as.data.frame(colData(from)), metadata = data.frame(), pkeys=vector())
# })


#' Subset with brackets
#'
#' Subset with brackets
#'
#' @param x dreamletProcessedData object
#' @param i indeces to extract
#'
#' @rdname extract-methods
#' @aliases [,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("[", signature(x="dreamletProcessedData"),
	function(x, i){   
		res = new("dreamletProcessedData", x@.Data[i], 
			data = x@data, 
			metadata = x@metadata,
			pkeys = x@pkeys)
		names(res) = names(x)[i]
		res
	}
)


setGeneric('assayNames', SummarizedExperiment::assayNames)
setGeneric('assay', SummarizedExperiment::assay)
setGeneric('colData', SummarizedExperiment::colData)
setGeneric('metadata', S4Vectors::metadata)

#' Get assayNames
#' 
#' Get assayNames
#' 
#' @param x dreamletProcessedData object
#' @param ... other arguments
#'
#' @rdname assayNames-methods
#' @aliases assayNames,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("assayNames", signature(x="dreamletProcessedData"),
	function(x, ...){   
		names(x)
	}
)

#' Get assay
#' 
#' Get assay
#' 
#' @param x dreamletProcessedData object
#' @param i number indicating index, or string indicating assay
#' @param withDimnames not used
#' @param ... other arguments
#'
#' @rdname assay-methods
#' @aliases assay,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("assay", signature(x="dreamletProcessedData"),
	function(x, i, withDimnames=TRUE,...){   
		x[[i]]
	}
)


#' Extract colData from dreamletProcessedData
#' 
#' Extract colData from dreamletProcessedData
#'
#' @param x A dreamletProcessedData object
#' @param ... other arguments
#' @export
setMethod("colData", "dreamletProcessedData",
	function(x,...){
		x@data
})


#' Extract metadata from dreamletProcessedData
#' 
#' Extract metadata from dreamletProcessedData
#'
#' @param x A dreamletProcessedData object
#' @export
setMethod("metadata", "dreamletProcessedData",
	function(x){
		x@metadata
})


#' Show object
#' 
#' Show object
#' 
#' @param object dreamletProcessedData object
#'
#' @rdname show-methods
#' @aliases show,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("show", "dreamletProcessedData",
	function(object){
		print(object)
	}
)



#' Print object
#' 
#' Print object
#' 
#' @param x dreamletProcessedData object
#' @param ... other arguments
#' 
#' @importFrom utils head tail
#' @importFrom S4Vectors coolcat
#' @export
#' @rdname print-methods
#' @aliases print,dreamletProcessedData,dreamletProcessedData-method
setMethod("print", "dreamletProcessedData",
	function(x,...){

		cat('class:', class(x), '\n')

		# assay
	    nms <- assayNames(x)
	    if (is.null(nms))
	        nms <- character(length(assays(x, withDimnames=FALSE)))
	    coolcat("assays(%d): %s\n", nms)

		# colData
	    nms <- names(colData(x))
	    if (is.null(nms))
	        nms <- character(length(colData(x, withDimnames=FALSE)))
	    coolcat("colData(%d): %s\n", nms)

	    # metadata
	    nms <- names(metadata(x))
	    if (is.null(nms))
	        nms <- character(length(metadata(x, withDimnames=FALSE)))
	    coolcat("metadata(%d): %s\n", nms)

		df_count = lapply(x, function(obj) dim(obj))
		df_count = do.call(rbind, df_count)

		cat('Samples:\n min:', min(df_count[,2]), '\n max:', max(df_count[,2]))
		cat('\nGenes:\n min:', min(df_count[,1]), '\n max:', max(df_count[,1]), '\n')

		# metadata
	    nms <- names(details(x))
	    if (is.null(nms))
	        nms <- character(length(metadata(x, withDimnames=FALSE)))
	    coolcat("details(%d): %s\n", nms)
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
#' Extract details from dreamletProcessedData
#'
#' @param object A dreamletProcessedData object
#'
#' @rdname details-methods
#' @export
setGeneric('details', function(object){
	standardGeneric("details")
	})


#' @export
#' @rdname details-methods
#' @aliases details,dreamletProcessedData-method
setMethod("details", "dreamletProcessedData",
	function(object){
			
	df = lapply( assayNames(object), function(k){

		obj = assay(object, k)
		DataFrame(assay = k, n_retained = ncol(obj$E), formula = paste(as.character(obj$formula), collapse=''))
	})
	df = do.call(rbind, df)

	df
})




#' @export
#' @rdname details-methods
#' @aliases details,dreamletResult-method
setMethod("details", "dreamletResult",
	function(object){
					
		object@df_details
})













