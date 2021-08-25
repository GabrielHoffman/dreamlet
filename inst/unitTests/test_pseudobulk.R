


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


# Check dreamlet::aggregateToPseudoBulk()
# compared to muscat::aggregateData
test_aggregateData = function(){

	library(muscat)

	# pseudobulk counts by cluster-sample
	data(example_sce)
	pb <- aggregateData(example_sce)

	library(SingleCellExperiment)
	assayNames(example_sce)  # one sheet per cluster
	head(assay(example_sce)) # n_genes x n_samples

	# scaled CPM
	cpm <- edgeR::cpm(assay(example_sce))
	assays(example_sce)$cpm <- cpm
	pb <- aggregateData(example_sce, assay = "cpm", scale = TRUE)
	# head(assay(pb)) 

	pb2 <- dreamlet::aggregateToPseudoBulk(example_sce, assay = "cpm", scale = TRUE,BPPARAM = SnowParam(2, progressbar=TRUE))
	
	check1 = checkEquals(pb, pb2)


	# aggregate by cluster only
	pb <- muscat::aggregateData(example_sce, by = "cluster_id")
	length(assays(pb)) # single assay
	head(assay(pb))    # n_genes x n_clusters

	pb2 <- dreamlet::aggregateToPseudoBulk(example_sce, by = "cluster_id")
	

	check1 & checkEquals(pb, pb2)
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


