# Gabriel Hoffman
# Jan 18, 2022


# Convert DelayedMatrix to sparseMatrix
# as(x, "sparseMatrix")
# da_to_sparseMatrix = function(x, verbose=FALSE){

# 	grid = colAutoGrid(x)

# 	res = lapply( seq(1,ncol(grid)), function(j, x, grid){
# 		if(verbose) cat("\r", j, " / ", ncol(grid), '    ')
# 		viewport <- grid[[1L, as.integer(j)]]

# 		.read_matrix_block(x, viewport)
# 		}, x=x, grid=grid)

# 	cat("\n")
# 	spMat <- cbind_list_of_sparseMatrix(res, verbose)
# 	rownames(spMat) = rownames(x)
# 	colnames(spMat) = colnames(x)

# 	spMat
# }
