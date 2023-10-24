
test_cellCounts = function(){

	data(example_sce)

	# create pseudobulk for each sample and cell cluster
	pb <- aggregateToPseudoBulk(example_sce, 
	   assay = "counts",    
	   cluster_id = 'cluster_id', 
	   sample_id = 'sample_id',
	   verbose=FALSE)

	counts1 = cellCounts(pb)

	counts2 = computeCellCounts(example_sce, 'cluster_id', 'sample_id')

	checkEquals(counts1, counts2)
}