# Gabriel Hoffman
# April 6, 2021



#' Remove constant terms from formula
#' 
#' Remove constant terms from formula
#'
#' @param formula original formula
#' @param data data.frame
#'
#' @examples
#' 
#' # Valid formula
#' removeConstantTerms(~group + extra, sleep)
#' 
#' # there is no variation in 'group' in this dataset
#' removeConstantTerms(~group + extra, sleep[1:3,])
#' 
#' @details Adapted from \code{MoEClust::drop_constants}
#' 
#' @importFrom stats terms update.formula reformulate as.formula
#' @importFrom lme4 findbars
#' @export
removeConstantTerms = function(formula, data){

	stopifnot(is(formula, "formula"))
	stopifnot(is(data, "data.frame"))

	# throw error if variable is not in data
	checkFormula( formula, data)

	trmf = terms(formula)

	intercept = attr(trmf, "intercept")
	fterms = attr(trmf, "term.labels")

	# find terms with bars, and add back parentheses
	i = which(fterms %in% findbars(formula)[])
	fterms[i] = paste0('(', fterms[i], ')')

	# keep only columns corresponding to variables in formula
	idx = colnames(data) %in% all.vars(update.formula(formula, NULL ~ .))
	data = data[,idx, drop=FALSE]

	# identify variables with no variation
	excludeVar = names(which(apply(data, 2L, function(x) all(x == x[1L], na.rm=TRUE))))
	
	# replace each excluded variable with an intercept
	if( length(excludeVar) > 0){
		fterms_new = array(sapply(excludeVar, function(x) gsub(x, '1', fterms)))
		fterms_new = fterms_new[fterms_new!=""]
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


# data = data.frame(df[1:4,], group_id=pb$group_id)[1:2,]

# formula = Size ~ (1|group_id) + Size +0

# removeConstantTerms(formula, data)












