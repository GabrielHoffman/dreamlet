# Gabriel Hoffman
# March 22, 2023
#
# Compute normalized counts

#' Compute normalized counts
#' 
#' Compute normalized counts as counts per million
#' 
#' @param sce \code{SingleCellExperiment} with counts stored as \code{counts(sce)}
#' 
#' @details This function gives same result as \code{edgeR::cpm(counts(sce), log=FALSE)}
#' 
#' @return matrix of CPM values
#' @seealso also \code{edgeR::cpm()}
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#' 
#' data(example_sce)
#' 
#' normcounts(example_sce) = computeNormCounts(example_sce)
#' @importFrom DelayedMatrixStats colSums2
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix t
#' @export
computeNormCounts = function(sce){
	# Compute library size
	lib.size = colSums2(counts(sce)) 

	t(t(counts(sce))/(lib.size/ 1e6))
}

#' Compute log normalized counts
#' 
#' Compute normalized counts as log2 counts per million
#' 
#' @param sce \code{SingleCellExperiment} with counts stored as \code{counts(sce)}
#' @param prior.count average count to be added to each observation to avoid taking log of zero
#' 
#' @details This function gives same result as \code{edgeR::cpm(counts(sce), log=TRUE)}
#' 
#' @return matrix of log CPM values
#' @seealso also \code{edgeR::cpm()}
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#' 
#' data(example_sce)
#' 
#' logcounts(example_sce) = computeLogCPM(example_sce)
#' @importFrom DelayedMatrixStats colSums2
#' @importFrom SingleCellExperiment counts
#' @export
computeLogCPM = function(sce, prior.count = 2){

	# Compute library size
	lib.size = colSums2(counts(sce)) 

	# compute prior count scaled by library size 
	# as in edgeR::cpm() call to C++ function add_prior.cpp
	pc = prior.count * lib.size / mean(lib.size)
	
	t(log2(t(counts(sce)) + pc) - log2(lib.size) + log2(1e6))
}