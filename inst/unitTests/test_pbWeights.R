




# test_pbWeights = function(){

# 	library(muscat)
# 	library(dreamlet)

# 	data(example_sce)

# 	# create pseudobulk for each sample and cell cluster
# 	pb <- aggregateToPseudoBulk(example_sce,
# 	  assay = "counts",
# 	  sample_id = "sample_id",
# 	  cluster_id = "cluster_id",
# 	  verbose = FALSE
# 	)

# 	# cell count weights
# 	####################
# 	weightsList <- pbWeights(example_sce,
# 	    sample_id = "sample_id",
# 	    cluster_id = "cluster_id")

# 	# use externallly comptuted cell weights
# 	res.proc1 <- processAssays(pb, ~ group_id, weightsList = weightsList)

# 	# use internally comptuted cell weights
# 	res.proc2 <- processAssays(pb, ~ group_id)

# 	for( i in seq(length(weightsList)) ){
# 		weightsList[[i]][] = 1
# 	}
# 	res.proc3 <- processAssays(pb, ~ group_id, weightsList = weightsList)

# 	# These should be identical, since weights of 1 are used
# 	checkIdentical(res.proc3, res.proc2)

# 	# delta weights
# 	###############
# 	weightsList <- pbWeights(example_sce,
# 	    sample_id = "sample_id",
# 	    cluster_id = "cluster_id", 
# 	    method = "delta")

# 	# use externallly comptuted cell weights
# 	res.proc1 <- processAssays(pb, ~ group_id, weightsList = weightsList)

# 	# use internally comptuted cell weights
# 	res.proc2 <- processAssays(pb, ~ group_id)

# 	# These should NOOOT be identical
# 	res = any(assay(res.proc1, 1)$weights == assay(res.proc2, 1)$weights)

# 	checkIdentical(res, FALSE)
# }