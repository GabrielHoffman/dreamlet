# Gabriel Hoffman
# April 6, 2021


#' Remove constant terms from formula
#' 
#' Remove constant terms from formula.  Also remove categorical variables with a max of one example per category
#'
#' @param formula original formula
#' @param data data.frame
#'
#' @return a formula, possibly with terms omitted.
#'
#' @examples
#' 
#' # Valid formula
#' removeConstantTerms(~ group + extra, sleep)
#' 
#' # there is no variation in 'group' in this dataset
#' removeConstantTerms(~ group + extra, sleep[1:3,])
#' 
#' @details Adapted from \code{MoEClust::drop_constants}
#' 
#' @importFrom stats terms update.formula reformulate as.formula
#' @importFrom lme4 findbars
#' @export
removeConstantTerms = function(formula, data){

	stopifnot(is(formula, "formula"))
	stopifnot(is(data, "data.frame"))

	# only retain columns used in the formula
	# Therefore, NA values in variable not used in the formula are ok
	data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]
	data = droplevels(data)

	# throw error if variable is not in data
	checkFormula( formula, data)

	trmf = terms(formula)

	intercept = attr(trmf, "intercept")
	fterms = attr(trmf, "term.labels")

	# if there are no terms, return the original formula
	if( length(fterms) == 0){
		return(formula)
	}

	# find terms with bars, and add back parentheses
	i = which(fterms %in% findbars(formula)[])
	fterms[i] = paste0('(', fterms[i], ')')

	# keep only columns corresponding to variables in formula
	idx = colnames(data) %in% all.vars(update.formula(formula, NULL ~ .))
	data = droplevels(data[,idx, drop=FALSE])

	# identify variables with no variation
	excludeVarConstant = names(which(apply(data, 2L, function(x){
		x = x[!is.na(x)]
		all(x == x[1L], na.rm=TRUE)
		})))
	
	# identify categorical variables with only single examples per category
	excludeVarCat = vapply(data, function(x){
		# exclude variable if it is a factor with max level count of 1
		# of no levels are observed more than once
		ifelse( is.factor(x) | is.character(x), (length(table(x)) == 1) | max(table(x)) == 1, FALSE)
		}, FUN.VALUE=logical(1))

	# combine excludes from multiple tests
	excludeVar = c(excludeVarConstant, names(excludeVarCat)[which(excludeVarCat)])
	excludeVar = unique(excludeVar)

	# replace each excluded variable with an intercept
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












