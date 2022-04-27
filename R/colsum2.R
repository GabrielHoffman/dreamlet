# Gabriel Hoffman
# Feb 22, 2022
#
# Faster version of DelayedArray::colsum

#' @import DelayedArray HDF5Array
colsum2 = function (x, group, reorder = TRUE, BPPARAM = SerialParam(), verbose=FALSE,...){

	# set block size
	# progressbar updates after each row chunk
	# rblock = ceiling(nrow(x)/10)
	# cblock = getAutoBlockSize() / (rblock*8)
	# grid = RegularArrayGrid(dim(x), spacings=c(rblock, cblock))

	.local <- function (x, group, reorder = TRUE, na.rm = FALSE, 
		grid = NULL,...){
		ugroup <- as.character(compute_ugroup(group, ncol(x), 
			reorder))
		if ( na.rm ) stop("na.rm = TRUE is not currently supported")
		# grid <- DelayedArray:::normarg_grid(grid, x)
		grid = defaultAutoGrid(x)
		block_results <- bplapply(seq_len(nrow(grid)), function(i) {
			.compute_colsum_for_grid_row(x, grid, i, group, ugroup, 
				na.rm = na.rm, verbose = verbose)
		}, BPPARAM=BPPARAM )
		ans <- do.call(rbind, block_results)
		dimnames(ans) <- list(rownames(x), ugroup)
		ans[,sort(colnames(ans)),drop=FALSE]
	}
	.local(x, group, reorder, ...)
}


#' @importFrom S4Vectors wmsg
compute_ugroup <- function(group, expected_group_len, reorder=TRUE){
	if (!(is.vector(group) || is.factor(group)))
		stop(wmsg("'group' must be a vector or factor"))
	if (length(group) != expected_group_len)
		stop(wmsg("incorrect length for 'group'"))
	if (!is.logical(reorder))
		stop(wmsg("'reorder' must be TRUE or FALSE"))
	## Taken from base::rowsum.default().
	ugroup <- unique(group)
	if (anyNA(ugroup))
		warning(wmsg("missing values for 'group'"))
	if (reorder)
		ugroup <- sort(ugroup, na.last=TRUE, method="quick")
	ugroup
}


.compute_colsum_for_grid_row <- function(x, grid, i, group, ugroup,
										 na.rm=FALSE, verbose=FALSE){
	grid_nrow <- nrow(grid)
	grid_ncol <- ncol(grid)
	ans <- matrix(0L, nrow=nrow(grid[[i, 1L]]), ncol=length(ugroup))
	## Inner loop on the grid cols. Sequential.
	for (j in seq_len(grid_ncol)) {
		if (verbose)
			message("\rProcessing block [[", i, "/", grid_nrow, ", ",
										   j, "/", grid_ncol, "]] ... ",
					appendLF=FALSE)
		block_ans <- .compute_colsum_for_block(x, grid, i, j,
											   group, na.rm=na.rm)
		m <- match(colnames(block_ans), ugroup)
		ans[ , m] <- ans[ , m] + block_ans
		if (verbose)
			message("OK", appendLF=FALSE)
	}
	ans
}

#' @importFrom S4Vectors extractROWS 
#' @importFrom IRanges ranges
.compute_colsum_for_block <- function(x, grid, i, j, group, na.rm=FALSE){
	viewport <- grid[[i, j]]
	block <- .read_matrix_block(x, viewport)
	group2 <- extractROWS(group, ranges(viewport)[2L])
	# colsum(block, group2, reorder=FALSE, na.rm=na.rm)
	colsum_beachmat(block, group2)
}




colsum_beachmat = function(x, group){

	if( is.factor(group) ){
		group = droplevels(group)
	}else if( is.list(group) ){

		grp = rep(0, max(vapply(group, max, numeric(1))))
		for( key in names(group) ){
			grp[group[[key]]] = key
		}
		group = factor(grp, sort(unique(grp)))

	}else{
		group = factor(group, sort(unique(group)))
	}

	if( is(x, "sparseMatrix") ){
		res = colsum_beachmat_sparseMatrix(x, group, unique(group))
	}else{
		res = colsum_beachmat_matrix(x, group, unique(group))
	}

	rownames(res) = rownames(x)
	colnames(res) = levels(group)

	res
}

