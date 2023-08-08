# Gabriel Hoffman
# Jan 10, 2021

create_idxlist <- function(fct) {
  fct <- droplevels(fct)

  res <- lapply(levels(fct), function(key) {
    which(fct == key)
  })
  names(res) <- levels(fct)
  res
}

# copied from DelayedArray:::.read_matrix_block
#' @importFrom methods as
#' @importFrom DelayedArray read_block
.read_matrix_block <- function(...) {
  block <- read_block(..., as.sparse = NA)
  if (is(block, "SparseArraySeed")) {
    block <- as(block, "CsparseMatrix")
  }
  block
}
