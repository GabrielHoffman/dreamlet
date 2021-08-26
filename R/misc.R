

# extract table of cell counts from 'int_colData'
# of pseudobulks as returned by 'aggregateData'
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment int_colData
.n_cells <- function(x) {
    y <- int_colData(x)$n_cells
    if (is.null(y)) return(NULL)
    if (length(metadata(x)$agg_pars$by) == 2)
        y <- as.matrix(data.frame(y, check.names = FALSE))
    return(as.table(y))
}


#' Check variables in a formula
#'
#' Check that variables in formula are present in the data
#'
#' @param formula formula of variables to check
#' @param data data.frame storing variables in the formula
#'
#' @importFrom lme4 subbars
#' @importFrom stats terms
checkFormula = function(formula, data){

	stopifnot(is(formula, "formula"))
	stopifnot(is(data, "data.frame"))

	v = attr(terms(subbars(formula)), "term.labels")
 	found = v %in% colnames(data)

 	if( any(!found) ){
 		stop("Variables in formula are not found in data:\n   ", paste(v[!found], collapse=", "))
 	}
 }


















