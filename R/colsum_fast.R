# Gabriel Hoffman
# Jan 10, 2021

create_idxlist = function(fct){

	fct = droplevels(fct)

	res = lapply(levels(fct), function(key){
		which(fct == key)
		})
	names(res) = levels(fct)
	res
}

# copied from DelayedArray:::.read_matrix_block
#' @importFrom methods as
.read_matrix_block = function(...){
    block <- read_block(..., as.sparse = NA)
    if (is(block, "SparseArraySeed")) 
        block <- as(block, "CsparseMatrix")
    block
}

# # DelayedArray::colsum can be slow due to data manipulation in R
# #  Uses a lot of memory when `group` has many levels
# # Here, colsum_fast() reads DelayedMatrix into memory 
# # 	but does manipulation and processing in Rcpp 
# #   Data is stored and return as a sparseMatrix to reduce memory usage

# #' @import DelayedArray HDF5Array
# #' @importFrom BiocParallel bplapply
# #' @importFrom methods as
# colsum_fast = function(x, group, grid=NULL, BPPARAM=SerialParam()){

# 	if( ! is.factor(group) ){
# 		group = factor(group, sort(unique(group)))
# 	}

# 	# get 2D grid of the DelayedArray
# 	# grid <- DelayedArray:::normarg_grid(grid, x)
# 	if( is.null(grid) ){
# 		grid = defaultAutoGrid(x)
# 	}	

# 	if( BPPARAM$progressbar ){
# 		message("Processing each block...")
# 	}

# 	# loop over chunks of columns
# 	result = bplapply( seq(1,ncol(grid)),function(j,x, grid, group){
# 		#cat("\r", j, '/', ncol(grid),'   ')

# 		# load within bplapply 
# 		suppressPackageStartupMessages(requireNamespace("DelayedArray"))
# 		suppressPackageStartupMessages(requireNamespace("HDF5Array"))

# 		# get indecies of column chunk to extract
# 		vp = grid[[1L, as.integer(j)]]
# 		# vp = DelayedArray:::getArrayElement(grid, c(1L, as.integer(j)))
# 		idx1 = start(vp@ranges)[2]
# 		idx2 = end(vp@ranges)[2]

# 		# get list of column assignments based on group
# 		idxlst = create_idxlist(group[seq(idx1, idx2)])

# 		# loop over chunks of row
# 		res = lapply( seq(1,nrow(grid)),function(i, x, grid){

# 			# get chunk of data
# 			viewport <- grid[[as.integer(i), as.integer(j)]]
# 			# viewport = DelayedArray:::getArrayElement(grid, c(1L, as.integer(j)))
# 			block = .read_matrix_block(x, viewport)

# 			# oringal (slow) code
# 			# res2 = DelayedArray::colsum(block, group[seq(idx1, idx2)])

# 			# Use RcppEigen to summarize data
# 			if( is(block, "sparseMatrix") ){	
# 				res3 = rowSums_by_chunk_sparse(block, idxlst, FALSE)
# 			}else{
# 				res3 = rowSums_by_chunk(block, idxlst, FALSE)
# 			}
# 			colnames(res3) = names(idxlst)

# 			# ensure that all chunks have same column sorting
# 			res3[,sort(colnames(res3)),drop=FALSE]
# 		}, x=x, grid=grid)
# 		if( length(res) > 1){
# 			# Aggregate row chunks
# 			res = do.call(rbind, res)
# 		}else{
# 			res = res[[1]]
# 		}
# 		# convert to sparseMatrix
# 		res = as(res, "sparseMatrix")

# 		# create fully expanded sparseMatrix
# 		# In the future, might be worth doing this trasnform inside 
# 		#    rowSums_by_chunk_sparse() to avoid passing data back to R
# 		# spMat = aggregateByColnames(list(res), list(colnames(res)), levels(group))
# 		spMat = aggregateByColnames1(res, colnames(res), levels(group))
# 		rownames(spMat) = rownames(x)
# 		colnames(spMat) = levels(group)

# 		spMat

# 	}, x=x, grid=grid, group=group, BPPARAM=BPPARAM)

# 	if( BPPARAM$progressbar ){
# 		message("Summing results from each block...")
# 	}

# 	# sum list of sparse matricies
# 	matSum = sumSpMatList(result, BPPARAM$progressbar)
# 	rownames(matSum) = rownames(result[[1]])
# 	colnames(matSum) = colnames(result[[1]])

# 	matSum
# }







	# Aggregate colums chunks by summing over column chunks
	# aggrResult = sparseMatrix(i=1, j=1, x=0, dims=c(nrow(x), nlevels(group)))
	# colnames(aggrResult) = levels(group)
	# rownames(aggrResult) = rownames(x)

	# for(batchResult in result){

	# 	idx = match(colnames(batchResult), colnames(aggrResult))

	# 	aggrResult[,idx] = aggrResult[,idx] + batchResult
	# }

	# aggrResult








