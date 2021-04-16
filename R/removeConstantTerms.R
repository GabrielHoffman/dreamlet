# Gabriel Hoffman
# April 6, 2021


#' Remove constant terms from formula
#' 
#' Remove constant terms from formula
#'
#' @param formula original formula
#' @param data data.frame
#' 
#' @importFrom lme4 subbars
#' @importFrom stats model.matrix as.formula update
#' @export
removeConstantTerms = function( formula, data){

	stopifnot(is(formula, "formula"))
	stopifnot(is(data, "data.frame"))

	# create design matrix
	design = tryCatch( {			
		model.matrix( subbars(formula), droplevels(data))
			}, 
		error = function(e) NULL)

	if( is.null(design) ) return(design)

	# which columns have zero variance
	isZero = apply(design, 2, var) == 0

	# identify assign to keep
	assignKeep = unique(c(0, attr(design,"assign")[which(!isZero)]))

	# find index of terms to exclude
	exclude = which(!attr(design,"assign") %in% assignKeep)

	if( length(exclude) == 0){
		return(formula)
	}

	# create part of formula to remove
	form_exclude = as.formula(paste(' ~ . -', paste(colnames(design)[exclude], collapse = ' - ')))

	# update formula to drop exclude terms
	update(formula, form_exclude)
}

