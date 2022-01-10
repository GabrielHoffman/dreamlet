






# create_idxlist = function(fct){

# 	fct = droplevels(fct)

# 	res = lapply(levels(fct), function(key){
# 		which(fct == key)
# 		})
# 	names(res) = levels(fct)
# 	res
# }

# # copied from DelayedArray:::.read_matrix_block
# .read_matrix_block = function(...){
#     block <- read_block(..., as.sparse = NA)
#     if (is(block, "SparseArraySeed")) 
#         block <- as(block, "CsparseMatrix")
#     block
# }



# #' @import BiocParallel DelayedArray
# fast_pb = function(x, group, grid=NULL, BPPARAM=SerialParam()){

# 	if( ! is.factor(group) ){
# 		group = factor(group, unique(group))
# 	}

# 	# get 2D grid of the DelayedArray
# 	grid <- DelayedArray:::normarg_grid(grid, x)

# 	# loop over chunks of columns
# 	result = bplapply( seq(1,ncol(grid)),function(j,x, grid, group){

# 		# load within bplapply 
# 	    suppressPackageStartupMessages(library("DelayedArray"))

# 	    # get indecies of column chunk to extract
# 		vp = grid[[1L, as.integer(j)]]
# 		# vp = DelayedArray:::getArrayElement(grid, c(1L, as.integer(j)))
# 		idx1 = start(vp@ranges)[2]
# 		idx2 = end(vp@ranges)[2]

# 		# get list of column assignments based on group
# 		idxlst = create_idxlist(group[seq(idx1, idx2)])

# 		# loop over chunks of row
# 		res = lapply( seq(1,nrow(grid)),function(i, x, grid){

# 			# get chunk of data
# 	    	viewport <- grid[[as.integer(i), as.integer(j)]]
# 			# viewport = DelayedArray:::getArrayElement(grid, c(1L, as.integer(j)))
# 			block = .read_matrix_block(x, viewport)

# 			# oringal (slow) code
# 			# res2 = DelayedArray::colsum(block, group[seq(idx1, idx2)])

# 			# Use RcppEigen to summarize data
# 			if( is(block, "sparseMatrix") ){	
# 				res3 = dreamlet:::rowSums_by_chunk_sparse(block, idxlst, FALSE)
# 			}else{
# 				res3 = dreamlet:::rowSums_by_chunk(block, idxlst, FALSE)
# 			}
# 			colnames(res3) = names(idxlst)

# 			# ensure that all chunks have same column sorting
# 			res3[,sort(colnames(res3))]
# 		}, x=x, grid=grid)
# 		# Aggregate row chunks
# 		res = do.call(rbind, res)
# 		# convert to sparseMatrix
# 		as(res, "sparseMatrix")
# 	}, x=x, grid=grid, group=group, BPPARAM=BPPARAM)


# 	# Aggregate colums chunks by summing over column chunks
# 	# aggrResult = sparseMatrix(i=1, j=1, x=0, dims=c(nrow(x), nlevels(group)))
# 	# colnames(aggrResult) = levels(group)
# 	# rownames(aggrResult) = rownames(x)

# 	# for(batchResult in result){

# 	# 	idx = match(colnames(batchResult), colnames(aggrResult))

# 	# 	aggrResult[,idx] = aggrResult[,idx] + batchResult
# 	# }

# 	# aggrResult
# }















