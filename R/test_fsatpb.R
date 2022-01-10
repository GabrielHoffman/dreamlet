






create_idxlist = function(a){

	unique_a = unique(a)
	res = lapply(unique_a, function(key){
		which(a == key)
		})
	names(res) = unique_a
	res
}

# copied from DelayedArray:::.read_matrix_block
.read_matrix_block = function(...){
    block <- read_block(..., as.sparse = NA)
    if (is(block, "SparseArraySeed")) 
        block <- as(block, "CsparseMatrix")
    block
}



#' @import BiocParallel DelayedArray
fast_pb = function(x, group, BPPARAM=SerialParam()){

	# get 2D grid of the DelayedArray
	grid <- DelayedArray:::normarg_grid(NULL, x)

	# loop over chunks of columns
	result = bplapply( seq(1,ncol(grid)),function(j,x, grid, group){

		# load within bplapply 
	    suppressPackageStartupMessages(library("DelayedArray"))

	    # get indecies of column chunk to extract
		vp = grid[[1L, as.integer(j)]]
		# vp = DelayedArray:::getArrayElement(grid, c(1L, as.integer(j)))
		idx1 = start(vp@ranges)[2]
		idx2 = end(vp@ranges)[2]

		# get list of column assignments based on group
		idxlst = create_idxlist(group[seq(idx1, idx2)])

		# loop over chunks of row
		res = lapply( seq(1,nrow(grid)),function(i, x, grid){

			# get chunk of data
	    	viewport <- grid[[as.integer(i), as.integer(j)]]
			# viewport = DelayedArray:::getArrayElement(grid, c(1L, as.integer(j)))
			block = .read_matrix_block(x, viewport)

			# oringal (slow) code
			# res2 = DelayedArray::colsum(block, group[seq(idx1, idx2)])

			# Use RcppEigen to summarize data
			if( is(block, "sparseMatrix") ){	
				res3 = rowSums_by_chunk_sparse(block, idxlst, FALSE)
			}else{
				res3 = rowSums_by_chunk(block, idxlst, FALSE)
			}
			colnames(res3) = names(idxlst)

			# ensure that all chunks have same column sorting
			res3[,sort(colnames(res3))]
		}, x=x, grid=grid)
		# Aggregate row chunks, and convert to sparseMatrix
		as(do.call(rbind, res), "sparseMatrix")
	}, x=x, grid=grid, group=group, BPPARAM=BPPARAM)

	# Aggregate colums chunks by summing
	Reduce('+', result )
}










