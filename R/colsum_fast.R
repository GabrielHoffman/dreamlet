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

# patched: https://github.com/GabrielHoffman/dreamlet/pull/23/commits/c453ac98ebc0329279b4dd3ae26a674df0e9b1f2
#' @importFrom methods as
#' @importFrom S4Arrays read_block
#' @importClassesFrom SparseArray SparseArray
.read_matrix_block <- function(...) {
  block <- read_block(..., as.sparse = NA)
  if (is(block, "SparseArray")) {
    block <- as(block, "CsparseMatrix")
  }
  block
}
