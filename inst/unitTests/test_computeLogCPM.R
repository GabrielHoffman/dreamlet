
test_computeNormCounts = function(){

	data(example_sce)

	# compute CPM using edgeR
	dge = DGEList(counts(example_sce))
	dge = calcNormFactors(dge, "none")

	value = max(abs(cpm(dge, log=FALSE) - computeNormCounts(example_sce) ))

	checkEqualsNumeric(value, 0)
}



test_computeLogCPM = function(){

	data(example_sce)

	# compute CPM using edgeR
	dge = DGEList(counts(example_sce))
	dge = calcNormFactors(dge, "none")

	prior.count = 1
	value = max(abs(cpm(dge, log=TRUE, prior.count=prior.count) - computeLogCPM(example_sce, prior.count=prior.count, scaledByLib=TRUE) ))

	checkEqualsNumeric(value, 0, tol=1e-2)
}


