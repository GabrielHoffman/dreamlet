
# TODO
# Estimateimate sigmaSq by ML
# fastest for DealyedArray

#' @export
getWeightFromCounts = function(countMatrix){

	count.gene = rowSums2(countMatrix) + 0.25
	count.lib = colSums2(countMatrix)
	ncell = ncol(countMatrix)
	sclSq = sum(count.lib^2)

	# var(countMatrix[j,] / count.lib)
	sigmaSq.hat.gene = rowVars(scale(countMatrix, scale = count.lib, center=FALSE), useNames=FALSE)
	sigmaSq.hat.gene[is.na(sigmaSq.hat.gene)] = 0

	# w = sapply( seq(nrow(countMatrix)), function(j){
	# 	1 / count.gene[j] + (sigmaSq.hat.gene[j] * sclSq) / (ncell^2 *count.gene[j]^2)
	# })

	# compute variance
	# vectorize
	v.hat = 1 / count.gene + (sigmaSq.hat.gene * sclSq) / (ncell^2 *count.gene^2)
	1 / v.hat
}

#' @export
getWeightsForCellType = function(sce, cluster_id, sample_id, CT){

	W = lapply( unique(sce[[sample_id]]), function(ID){

		idx = sce[[cluster_id]] == CT & sce[[sample_id]] == ID
		countMatrix = counts(sce)[,idx,drop=FALSE]

		getWeightFromCounts( countMatrix )
		})
	W = do.call(cbind, W)
	W = W / rowMeans2(W, useNames=FALSE)

	colnames(W) = unique(sce[[sample_id]])
	rownames(W) = rownames(sce)
	W
}

#' @export
getWeightsList = function(sce, cluster_id, sample_id){
	W.list = lapply( unique(sce[[cluster_id]]), function(CT){
		getWeightsForCellType( sce, cluster_id, sample_id, CT)
	})
	names(W.list) = unique(sce[[cluster_id]])
 	W.list
 }
