
#' Class dreamletProcessedData
#'
#' Class \code{dreamletProcessedData} 
#'
#' @name dreamletProcessedData-class
#' @rdname dreamletProcessedData-class
#' @exportClass dreamletProcessedData
#' @return none
setClass("dreamletProcessedData", contains="list", slots = c(data = 'data.frame', metadata='data.frame', by="vector"))

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
setMethod("[", signature(x="dreamletProcessedData"),
	function(x, i){   
		res = new("dreamletProcessedData", x@.Data[i], 
			data = x@data, 
			metadata = x@metadata, 
			by = x@by)
		names(res) = names(x)[i]
		res
	}
)


setGeneric('assayNames', SummarizedExperiment::assayNames)
setGeneric('assay', SummarizedExperiment::assay)
setGeneric('colData', SummarizedExperiment::colData)
setGeneric('colData<-', SummarizedExperiment::`colData<-`)
setGeneric('metadata', S4Vectors::metadata)

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
setMethod("assayNames", signature(x="dreamletProcessedData"),
	function(x, ...){   
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
setMethod("assay", signature(x="dreamletProcessedData"),
	function(x, i, withDimnames=TRUE,...){   
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
setMethod("colData", "dreamletProcessedData",
	function(x,...){
		x@data
})



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
setMethod("colData<-", "dreamletProcessedData",
    function(x, ..., value){
    
    # convert to data.frame
    value = as.data.frame(value)

    # check dimensions
    if( nrow(x@data) != nrow(value) ){
    	stop("Number of rows in colData(x) must remain the same")
    }

    # check same rownames
    if( ! all.equal(rownames(x@data), rownames(value)) ){    	
    	stop("rownames(colData(x)) must remain the same in the new data.frame")
    }

    x@data = value
    x
})

#' Extract metadata from \code{dreamletProcessedData}
#' 
#' Extract metadata from \code{dreamletProcessedData}
#'
#' @param x A dreamletProcessedData object
#'
#' @return object from \code{metadata} field
#' @aliases metadata,dreamletProcessedData,dreamletProcessedData-method
#' @export
setMethod("metadata", "dreamletProcessedData",
	function(x){
		x@metadata
})


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
setMethod("show", "dreamletProcessedData",
	function(object){

	cat('class:', class(object), '\n')

	# assay
    nms <- assayNames(object)
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    coolcat("assays(%d): %s\n", nms)

	# colData
    nms <- names(colData(object))
    if (is.null(nms))
        nms <- character(length(colData(object, withDimnames=FALSE)))
    coolcat("colData(%d): %s\n", nms)

    # metadata
    nms <- names(metadata(object))
    if (is.null(nms))
        nms <- character(length(metadata(object, withDimnames=FALSE)))
    coolcat("metadata(%d): %s\n", nms)

	df_count = lapply(object, function(obj) dim(obj))
	df_count = do.call(rbind, df_count)

	if( is.null(df_count) ){
		cat("No assays retained\n")
	}else{
		cat('Samples:\n min:', min(df_count[,2]), '\n max:', max(df_count[,2]))
		cat('\nGenes:\n min:', min(df_count[,1]), '\n max:', max(df_count[,1]), '\n')

		# metadata
	    nms <- names(details(object))
	    if (is.null(nms))
	        nms <- character(length(metadata(object, withDimnames=FALSE)))
	    coolcat("details(%d): %s\n", nms)
	}
})



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
setMethod("print", "dreamletProcessedData",
	function(x,...){
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
setGeneric('details', getGeneric("details", package="GSEABase"))


#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce, 
#'    assay = "counts",    
#'    cluster_id = 'cluster_id', 
#'    sample_id = 'sample_id',
#'    verbose=FALSE)
#'
#' # voom-style normalization
#' res.proc = processAssays( pb, ~ group_id)
#' 
#' # For each cell type, number of samples retained, 
#' # and variables retained
#' details(res.proc)
#' 
#' @export
#' @rdname details-methods
#' @aliases details,dreamletProcessedData-method
setMethod("details", "dreamletProcessedData",
	function(object){
			
	df = lapply( assayNames(object), function(k){

		obj = assay(object, k)

		if( is(obj, "EList") ){
			res = DataFrame(assay = k, n_retained = ncol(obj), formula = paste(as.character(obj$formula), collapse=''))
		}else{
			res = DataFrame(assay = k, n_retained = ncol(obj))
		}
		res
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





#' Extract expression and colData
#' 
#'
#' @param x A \code{dreamletProcessedData} object
#' @param assay assay to extract
#'
#' @rdname extractData-methods
#' @export
setGeneric('extractData', function(x, assay) standardGeneric("extractData"))


#' Extract expression and \code{colData}
#' 
#' Extract expression and \code{colData} from \code{dreamletProcessedData}
#'
#' @param x A \code{dreamletProcessedData} object
#' @param assay assay to extract
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
#'    assay = "counts",    
#'    cluster_id = 'cluster_id', 
#'    sample_id = 'sample_id',
#'    verbose=FALSE)
#'
#' # voom-style normalization
#' res.proc = processAssays( pb, ~ group_id)
#' 
#' # Extract data.frame of colData merged with expression.
#' # variables and genes are stored as columns, samples as rows
#' df_merge = extractData( res.proc, "B cells")
#' 
#' dim(df_merge)
#' 
#' # first few columns
#' df_merge[, 1:6]
#'
#' @importFrom S4Vectors merge
#' @rdname extractData-methods
#' @aliases extractData,dreamletProcessedData-method
#' @export
setMethod("extractData", c(x="dreamletProcessedData", assay="character"),
	function(x, assay){
	
	if( ! assay %in% assayNames(x) ){
		stop("assay not found: ", assay)
	} 				

	merge(colData(x), t(assay(x, assay)$E), by="row.names")
})










