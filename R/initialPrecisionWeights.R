
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
getWeightsForCellType = function(sce, cluster_id, sample_id, CT, weightCap){

	W <- lapply( unique(sce[[sample_id]]), function(ID){

		idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == ID
		countMatrix = counts(sce)[,idx,drop=FALSE]

		getWeightFromCounts( countMatrix )
		})
	W <- do.call(cbind, W)

	# scale each gene by min value
	# so lowest weight is now 1
	W <- W / rowMins( W, useNames=FALSE)

	# set a weight cap at weightCap
	W[W > weightCap] <- weightCap

	# scale each gene to have a mean of 1
	W <- W / rowMeans2(W, useNames=FALSE)

	colnames(W) <- unique(sce[[sample_id]])
	rownames(W) <- rownames(sce)
	W
}

#' @export
getWeightsList = function(sce, cluster_id, sample_id, weightCap = 10){

	if( ! cluster_id %in% colnames(colData(sce)) ){
		msg <- paste0("sample_id entry not found in colData(sce): ", cluster_id)
		stop( msg )
	}
	if( ! sample_id %in% colnames(colData(sce)) ){
		msg <- paste0("sample_id entry not found in colData(sce): ", sample_id)
		stop( msg )
	}

	W.list <- lapply( unique(sce[[cluster_id]]), function(CT){
		getWeightsForCellType( sce, cluster_id, sample_id, CT, weightCap)
	})
	names(W.list)<- unique(sce[[cluster_id]])
 	W.list
}



#' @export
trimWeightOutliersGene = function(x, zmax){

	# compute z-score
	zscore = scale(x)

	# extract parameters of transform
	# z-score = (x - mu) / s
	mu = attr(zscore,"scaled:center")
	s = attr(zscore,"scaled:scale")
	
	# if x exceeds original value giving z-score of zmax, 
	# replace with that orginal value
	x[x > zmax * s + mu] = zmax * s + mu

	# normalize values to have a mean of 1
	x / mean(x)
}

#' @export
trimWeightOutliers = function(W, zmax = 5){

	t(apply(W, 1, trimWeightOutliersGene, zmax = zmax))
}


