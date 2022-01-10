


library(RUnit)

# Check dreamlet:::summarizeAssayByGroup2()
# compared to scuttle::summarizeAssayByGroup()
test_pseudobulk_example = function(){
                
	example_sce <- scuttle::mockSCE()

	ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)

	out <- scuttle::summarizeAssayByGroup(example_sce, ids)

	out2 <- dreamlet:::summarizeAssayByGroup2(example_sce, ids)

	checkEquals(out, out2)
}

test_rowSums_by_chunk = function(){

	set.seed(17)# to be reproducible
	n = 400
	p = 1000
	M <- Matrix::rsparsematrix(n, p, density=.1)

	idxlist = list(1:p)

	res = dreamlet:::rowSums_by_chunk(as.matrix(M), idxlist, TRUE)
	res2 = dreamlet:::rowSums_by_chunk_sparse(M, idxlist, TRUE)

	checkEqualsNumeric(res, res2)
}



# Check dreamlet::aggregateToPseudoBulk()
# compared to muscat::aggregateData
test_aggregateData = function(){

	library(muscat)
	library(BiocParallel)

	# pseudobulk counts by cluster-sample
	data(example_sce)
	pb <- aggregateData(example_sce)

	library(SingleCellExperiment)
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
	
	check1 = checkEquals(pb, pb2)

	# aggregate by cluster only
	pb <- muscat::aggregateData(example_sce, by = "cluster_id")
	# length(assays(pb)) # single assay
	# head(assay(pb))    # n_genes x n_clusters

	pb2 <- dreamlet::aggregateToPseudoBulk(example_sce, cluster_id = "cluster_id")

	check1 & checkEquals(pb, pb2)
}

test_aggregateToPseudoBulk_datatype = function(){

	# compare pseudobulk by rowSums from DelayedMatrix, matrix, and sparseMatrix
	library(muscat)
	library(SingleCellExperiment)
	library(DelayedArray)

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
			sample_id = "sample_id") 

	checkEquals(pb_sparseMatrix, pb_matrix ) & checkEquals(pb_sparseMatrix, pb_delayedMatrix ) 
}


test_pmetadata = function(){


	# devtools::reload("/Users/gabrielhoffman/workspace/repos/dreamlet")

	library(muscat)
	library(BiocParallel)
	library(SingleCellExperiment)

	# pseudobulk counts by cluster-sample
	data(example_sce)
	pb <- aggregateData(example_sce[1:100,])

	# simulated pmetadata at the assay level
	df = data.frame( ID = colnames(pb), assay = sort(rep(assayNames(pb), ncol(pb))))
	df$Size = rnorm(nrow(df))
	pkeys = c("ID", "assay")

	res.proc = processAssays( pb, ~ (1|group_id) + Size, min.count=5, pmetadata=df, pkeys=pkeys)

	fig = plotVoom(res.proc)

	vp.lst = fitVarPart( res.proc, ~ Size + (1|group_id))

	fig = plotVarPart( sortCols(vp.lst))

	res.dl = dreamlet( res.proc, ~ (1|group_id) + Size)

	tab = topTable(res.dl, coef="Size")

	TRUE
}


# devtools::reload("/Users/gabrielhoffman/workspace/repos/dreamlet")

# pb2 <- dreamlet::aggregateToPseudoBulk(example_sce, assay = "cpm", scale = TRUE,BPPARAM = SnowParam(2, progressbar=TRUE))



# out2 <- dreamlet:::summarizeAssayByGroup2(example_sce, ids,BPPARAM = SnowParam(2, progressbar=TRUE))

# Loading required package: HDF5Array
# Loading required package: DelayedArray
# Loading required package: stats4
# Loading required package: Matrix
# Loading required package: BiocGenerics
# Loading required package: parallel
# rhdf5
# HDF5Array


