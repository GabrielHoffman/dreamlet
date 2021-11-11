

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
#' @return If formula is valid, return TRUE.  Else throw error
#'
#' @examples
#'
#' # Valid formula
#' dreamlet:::checkFormula( ~ speed, cars)
#' 
#' # Not valid formula
#' # dreamlet:::checkFormula( ~ speed + a, cars)
#' 
#' @importFrom stats terms
checkFormula = function(formula, data){

	stopifnot(is(formula, "formula"))
	stopifnot(is(data, "data.frame"))

	v = all.vars(formula)
 	found = v %in% colnames(data)

 	if( any(!found) ){
 		stop("Variables in formula are not found in data:\n   ", paste(v[!found], collapse=", "))
 	}
 }






#' Check if two formulas are equal
#'
#' Check if two formulas are equal by evaluating the formulas and extracting terms
#'
#' @param formula1 first formula
#' @param formula2 second formula
#'
#' @return boolean value indciating of formulas are equivalent
#'
#' @examples
#'
#' # These formulas are equivalent
#' formula1 = ~ Size + 1
#' formula2 = ~ 1 + Size 
#'
#' dreamlet:::equalFormulas( formula1, formula2)
#'
#' @importFrom stats terms
equalFormulas = function(formula1, formula2){

    # extract terms from forumula1
    trmf1 = terms(as.formula(formula1))
    intercept1 = attr(trmf1, "intercept")
    fterms1 = attr(trmf1, "term.labels")

    # extract terms from forumula2
    trmf2 = terms(as.formula(formula2))
    intercept2 = attr(trmf2, "intercept")
    fterms2 = attr(trmf2, "term.labels")

    # check equality of intercept and variables 
    # sort because the order might be different
    (intercept1 == intercept2) & identical(sort(fterms1), sort(fterms2))
}


#' Extract residuls from dreamletResult
#' 
#' Extract residuls from dreamletResult
#' 
#' @param object dreamletResult object
#' @param ... other arguments
#' 
#' @rdname residuals-methods
#' @aliases residuals,dreamletResult,dreamletResult-method
setMethod("residuals", "dreamletResult",
  function(object,...){

    res = lapply( names(object), function(key){
        residuals(object[[key]])
    })
    names(res) = names(object)

    res
})












