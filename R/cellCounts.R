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
#' @export
cellCounts = function(x){

	if( ! is(x, "SingleCellExperiment") ){
		stop("x must be SingleCellExperiment")
	}

	do.call(rbind, int_colData(x)$n_cells)
}



