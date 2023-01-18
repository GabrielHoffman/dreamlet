# Gabriel Hoffman
# Oct 28, 2021

#' Extract cell counts
#' 
#' Extract matrix of cell counts from \code{SingleCellExperiment}
#' 
#' @param x a \code{SingleCellExperiment}
#' 
#' @return matrix of cell counts with samples as rows and cell types as columns
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
#' # get matrix of cell counts for each sample
#' cellCounts(pb)
#'
#' @seealso \code{computeCellCounts()}
#' @export
cellCounts = function(x){

	if( ! is(x, "SingleCellExperiment") ){
		stop("x must be SingleCellExperiment")
	}

	do.call(rbind, int_colData(x)$n_cells)
}


#' Get cell counts with metadata
#' 
#' Get cell counts with metadata for each sample
#'  
#' @param sce \code{SingleCellExperiment}
#' @param annotation string indicating column in \code{colData(sce)} storing cell type annotations
#' @param sampleIDs string indicating column in \code{colData(sce)} storing sample identifers
#' 
#' @return \code{matrix} storing cell counts
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' counts = computeCellCounts(example_sce, "cluster_id", "sample_id")
#'
#' counts[1:4, 1:4]
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom dplyr `%>%` group_by across summarize distinct n as_tibble
#' @export
computeCellCounts = function(sce, annotation, sampleIDs){

	# count number of each cell type observed for each sample
	df = colData(sce) %>%
		as_tibble %>%
		group_by(across(annotation), across(sampleIDs)) %>%
		summarize( count=n(), .groups='keep') %>%
		distinct 

	if( ! is.factor(df[[annotation]])){
		df[[annotation]] = factor(df[[annotation]])
	}
	if( ! is.factor(df[[sampleIDs]])){
		df[[sampleIDs]] = factor(df[[sampleIDs]])
	}

	# convert to matrix form
	M = sparseMatrix( as.numeric(df[[annotation]]), 
	        as.numeric(df[[sampleIDs]]), 
	        x=df$count,
	        dims=c(nlevels(df[[annotation]]), nlevels(df[[sampleIDs]])))
	rownames(M) = levels(df[[annotation]])
	colnames(M) = levels(df[[sampleIDs]])

	t(as.matrix(M))
}















