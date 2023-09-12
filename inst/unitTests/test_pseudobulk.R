


# Check dreamlet:::summarizeAssayByGroup2()
# compared to scuttle::summarizeAssayByGroup()
test_pseudobulk_example = function(){
                
	example_sce <- scuttle::mockSCE()

	ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)

	out <- scuttle::summarizeAssayByGroup(example_sce, ids)

	out2 <- dreamlet:::summarizeAssayByGroup2(example_sce, ids, statistics = c("mean", "sum", "num.detected", "prop.detected", "median"))

	checkEquals(out, out2)
}


# Check dreamlet::aggregateToPseudoBulk()
# compared to muscat::aggregateData
test_aggregateData = function(){

	# pseudobulk counts by cluster-sample
	data(example_sce)
	pb <- aggregateData(example_sce)

	# assayNames(example_sce)  # one sheet per cluster
	# head(assay(example_sce)) # n_genes x n_samples

	# scaled CPM
	cpm <- edgeR::cpm(assay(example_sce))
	assays(example_sce)$cpm <- cpm
	pb <- aggregateData(example_sce, assay = "cpm", scale = TRUE)
	# head(assay(pb)) 

	pb2 <- dreamlet::aggregateToPseudoBulk(example_sce, assay = "cpm",
			cluster_id = "cluster_id",
			sample_id = "sample_id", 
			scale = TRUE) 

	metadata(pb2)$aggr_means = c()

	# convert from sparseMatrix to matrix just for comparing to expectation
	for(id in assayNames(pb2)){
		assay(pb2, id) = as.matrix(assay(pb2, id))
	}
	
	check1 = checkEquals(pb, pb2)

	# aggregate by cluster only
	pb <- muscat::aggregateData(example_sce, by = "cluster_id")
	# length(assays(pb)) # single assay
	# head(assay(pb))    # n_genes x n_clusters

	pb2 <- dreamlet::aggregateToPseudoBulk(example_sce, cluster_id = "cluster_id")
	assays(pb2)[[1]] = as.matrix(assays(pb2)[[1]])
	metadata(pb2)$aggr_means = c()

	check1 & checkEquals(pb, pb2)
}

test_colsum_fast = function(){

	set.seed(17)# to be reproducible
	n = 400
	p = 500
	M <- rsparsematrix(n, p, nnz = n*p*.1)
	# M[] = 1L

	group = as.character(sample.int(200, p, replace=TRUE))

	res1 = dreamlet:::colsum2(DelayedArray(M), group)

	res2 = DelayedArray::colsum(DelayedArray(M), group)

	RUnit::checkEquals(as.matrix(res1), res2 )
}

test_aggregateToPseudoBulk_datatype = function(){

	# compare pseudobulk by rowSums from DelayedMatrix, matrix, and sparseMatrix

	# pseudobulk counts by cluster-sample
	data(example_sce)

	is(assay(example_sce, "counts"))

	# sparseMatrix
	pb_sparseMatrix <- dreamlet::aggregateToPseudoBulk(example_sce, assay = "counts",
			cluster_id = "cluster_id",
			sample_id = "sample_id") 

	# matrix
	assay(example_sce, "counts") = as.matrix(assay(example_sce, "counts"))
	pb_matrix <- dreamlet::aggregateToPseudoBulk(example_sce, assay = "counts",
			cluster_id = "cluster_id",
			sample_id = "sample_id") 

	# matrix
	assay(example_sce, "counts") = DelayedArray(assay(example_sce, "counts"))
	pb_delayedMatrix <- dreamlet::aggregateToPseudoBulk(example_sce, assay = "counts",
			cluster_id = "cluster_id",
			sample_id = "sample_id",
			BPPARAM=SnowParam(2, progressbar=TRUE)) 

	# convert from sparseMatrix to matrix just for comparing to expectation
	# for(id in assayNames(pb_delayedMatrix)){
	# 	assay(pb_delayedMatrix, id) = as.matrix(assay(pb_delayedMatrix, id))
	# }

	assay(pb_sparseMatrix, 1)[1:3, 1:3]
	assay(pb_matrix, 1)[1:3, 1:3]
	assay(pb_delayedMatrix, 1)[1:3, 1:3]

	checkEquals(pb_sparseMatrix, pb_matrix ) & checkEquals(pb_sparseMatrix, pb_delayedMatrix ) 
}

test_cell_level_means = function(){

	data(example_sce)

	# continuous value
	example_sce$value1 = rnorm(ncol(example_sce))
	example_sce$value2 = rnorm(ncol(example_sce))

	# create pseudobulk for each sample and cell cluster
	pb <- aggregateToPseudoBulk(example_sce, 
	   assay = "counts",    
	   cluster_id = 'cluster_id', 
	   sample_id = 'sample_id',
	   verbose=FALSE)

	metadata(pb)$aggr_means

	metadata(pb)$agg_pars$by

	# voom-style normalization
	res.proc = processAssays( pb, ~ group_id + value1)

	vp = fitVarPart(res.proc, ~ group_id + value1)

	fit = dreamlet(res.proc, ~ group_id + value1)
}



# test_da_to_sparseMatrix 