
# TODO
# Estimateimate sigmaSq by ML
# fastest for DealyedArray

#' @export
getWeightFromCounts = function(countMatrix){

	count.gene <- rowSums2(countMatrix, useNames=FALSE)
	count.lib <- colSums2(countMatrix, useNames=FALSE)
	ncell <- ncol(countMatrix)
	sclSq <- sum(count.lib^2)

	# add pseudocount
	count.gene <- count.gene + 0.25
	count.lib <- count.lib + 1

	# normalize counts by library size
	# add pseudocount to counts here
	normCounts <- scale(countMatrix + 0.25, 
					scale = count.lib, 
					cente = FALSE)
	# compute variance for each row
	sigmaSq.hat.gene <- rowVars(normCounts, useNames=FALSE)
	sigmaSq.hat.gene[is.na(sigmaSq.hat.gene)] <- 0

	# compute variance
	# vectorize
	v.hat <- 1 / count.gene + (sigmaSq.hat.gene * sclSq) / (ncell^2 *count.gene^2)
	1 / v.hat
}

#' @export
getWeightsForCellType = function(sce, cluster_id, sample_id, CT){

	W <- lapply( unique(sce[[sample_id]]), function(ID){

		idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == ID
		countMatrix = counts(sce)[,idx,drop=FALSE]

		getWeightFromCounts( countMatrix )
		})
	W <- do.call(cbind, W)
	W <- W / rowMeans2(W, useNames=FALSE)

	colnames(W) <- unique(sce[[sample_id]])
	rownames(W) <- rownames(sce)
	W
}

#' @export
getWeightsList = function(sce, cluster_id, sample_id){
	W.list <- lapply( unique(sce[[cluster_id]]), function(CT){
		getWeightsForCellType( sce, cluster_id, sample_id, CT)
	})
	names(W.list)<- unique(sce[[cluster_id]])
 	W.list
 }
