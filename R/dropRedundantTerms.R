
#' Drop redundant terms from the model
#' 
#' Detect co-linear fixed effects and drop the last one
#'
#' @param formula original formula
#' @param data data.frame
#' @param tol tolerance to test difference of correlation from 1 or -1
#'
#' @return a formula, possibly with terms omitted.
#'
#' @examples
#' 
#' # Valid formula
#' dropRedundantTerms(~ group + extra, sleep)
#' 
#' @importFrom stats terms update.formula reformulate as.formula cor
#' @importFrom lme4 findbars
#' @importFrom Matrix summary
#' @export
dropRedundantTerms = function(formula, data, tol=1e-3){

	stopifnot(is(formula, "formula"))
	stopifnot(is(data, "data.frame"))

	# only retain columns used in the formula
	# Therefore, NA values in variable not used in the formula are ok
	data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]
	data = droplevels(data)

	# throw error if variable is not in data
	checkFormula( formula, data)

	# remove response
	form2 = update.formula(formula, NULL ~ .)

	# create design matrix
	dsgn = model.matrix(nobars(form2), data)

	# drop intercept
	idx = match("(Intercept)", colnames(dsgn))
	if( !is.na(idx)) dsgn = dsgn[,-idx,drop=FALSE]

	# identify redundant pairs of variables
	M = as(cor(dsgn), "sparseMatrix")
	M[is.na(M)] = 0
	M[lower.tri(M,diag=TRUE)] = 0
	df = summary(M)
	df = df[1-abs(df$x) < tol,,drop=FALSE]

	if( nrow(df) == 0){
		return(formula)
	}

	excludeVar = colnames(dsgn)[df$j]

	fterms = attr(terms(formula), "term.labels")

	# find terms with bars, and add back parentheses
	i = which(fterms %in% findbars(formula)[])
	fterms[i] = paste0('(', fterms[i], ')')

	if( length(excludeVar) > 0){

		# replace each exclude term with with intercept
		# gsub can overwrite other variables
		fterms_new = fterms
		for( x in excludeVar){

			# fterms_new = gsub(x, '1', fterms_new)

			# make sure string goes to end, or end followed by )
			fterms_new = gsub(paste0(x, '$'), '1', fterms_new)
			fterms_new = gsub(paste0(x, '\\)'), '1)', fterms_new)
		}

		fterms_new = unique(fterms_new)

		# fterms_new = array(sapply(excludeVar, function(x) gsub(x, '1', fterms)))
		fterms_new = fterms_new[fterms_new!=""]
		fterms_new = fterms_new[fterms_new!="1"]
	}else{
		fterms_new = fterms
	}

	# remove V:1 or 1:V or | 1
	exclude = grep(":1|1:|\\| 1", fterms_new)

	if( length(exclude) > 0){
		fterms_new = fterms_new[-exclude]
	}

	# get response to add back
	response = all.vars(update.formula(formula, . ~ NULL))
	if( attr(trmf, "response") == 0) response = NULL	

	# convert to formula
	if( length(fterms_new) > 0 ){
		form_new = reformulate(fterms_new, response=response)
	}else{
		form_new = ~ 1
	}

	# if no intercept in original formula, remove it here
	if(intercept == 0){
		form_new = update.formula(form_new, ~ . -1)       
	}

	environment(form_new) <- environment(formula)
	form_new
}

