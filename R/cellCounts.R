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
#'  
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
#' @export
cellCounts = function(x){

	if( ! is(x, "SingleCellExperiment") ){
		stop("x must be SingleCellExperiment")
	}

	do.call(rbind, int_colData(x)$n_cells)
}



