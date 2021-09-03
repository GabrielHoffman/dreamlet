

#' Test difference in cell type composition
#'
#' Test difference in cell type composition using a regression models
#'
#' @param obj \code{SingleCellExperiment} from \code{aggregateToPseudoBulk}
#' @param formula formula indicating variables to include in the regression
#' @param coef which coefficient to test
#' @param method which regression model to use
#'
#' @description Evaluates regression model where cell fraction is response using one of the 6 specified methods.  Methods can be fit with random effects, except for 'betabinomial'
#'
#' @export
cellTypeCompositionTest = function( obj, formula, coef, method = c("binomial", "betabinomial", "lm", "lmlog", "poisson", "nb") ){

	method = match.arg(method)

	# extract cell counts and other meta-data
	df_cellCount = do.call(rbind, int_colData(obj)$n_cells)
	df_cellCount = cbind(df_cellCount, TotalCells = rowSums(df_cellCount))
	df_cellCount = merge(df_cellCount, colData(obj), by="row.names")

  	# loop thru assays and perform test
	df_fit = lapply( assayNames(obj), function(k){

		# create formula
		if( method %in% c("binomial", "betabinomial") ){

			form = paste0("cbind(`", k, "`, TotalCells - `", k, "`) ~ ", as.character(formula)[2])

		}else if( method == "lm"){
			form = paste0("`", k, "`  / TotalCells ~ ", as.character(formula)[2])
		}else if(method == 'lmlog'){
			form = paste0("log(`", k, "`  / TotalCells) ~ ", as.character(formula)[2])
		}else if(method %in% c("poisson", "nb") ){		
			form = paste0("`", k, "` ~ offset(log(TotalCells)) + ", as.character(formula)[2])
		}

		# evaluate regression model
		res = regModel(form, df_cellCount, coef, method)

		data.frame(assay =k, res)
	})
	df_fit = as.data.frame(do.call(rbind, df_fit))
	colnames(df_fit)[3:5] = c("se", "zstat", "pValue")
	df_fit$p.adj = p.adjust(df_fit$pValue, "fdr")

	df_fit
}

.isMixedModelFormula = function(formula){	
    !is.null(findbars(as.formula(formula)))
}

#' @importFrom aod betabin
#' @importFrom MASS glm.nb
#' @importFrom lme4 glmer glmer.nb lmerControl glmerControl .makeCC
#' @importFrom lmerTest lmer
#' @importFrom stats lm glm glm.control
regModel = function(formula, data, coef, method = c("binomial", "betabinomial", "lm", "lmlog", "poisson", "nb")){

	method = match.arg(method)

	formula = as.formula(formula)

	if( .isMixedModelFormula(formula) ){

		# mixed model

		if( method %in% c("binomial",  "poisson") ){
			control = glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
			fit = glmer(formula, data, family=method, control=control)
			df = coef(summary(fit))[coef,,drop=FALSE]

		}else if(method == "betabinomial"){
			stop("'betabinomial' model can't accept random effects")

		}else if(method %in% c("lm", "lmlog") ){
			control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
	
			fit = lmer(formula, data, control = control)			
			colids = c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
			df = coef(summary(fit))[coef,colids,drop=FALSE]

		}else if(method == "nb" ){
			control = glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
			# if NB fit fails, use poisson model
			fit = tryCatch({				
				glmer.nb(formula, data, control=control)
				}, error = function(e){
				glmer(formula, data, family="poisson", control=control)
				}, warning=function(w) {				
				glmer(formula, data, family="poisson", control=control)
				})
			df = coef(summary(fit))[coef,,drop=FALSE]
		}
	}else{
		# only fixed effects
		if( method %in% c("binomial",  "poisson") ){
			fit = glm(formula, data, family=method)
			df = coef(summary(fit))[coef,,drop=FALSE]

		}else if(method == "betabinomial"){
			fit = betabin(formula, ~ 1, data)
			df = summary(fit)@Coef[coef,,drop=FALSE]

		}else if(method %in% c("lm", "lmlog") ){
			fit = lm(formula, data)
			df = coef(summary(fit))[coef,,drop=FALSE]

		}else if(method == "nb" ){
			# if NB fit fails, use poisson model
			fit = tryCatch({
				glm.nb(formula, data, control=glm.control(maxit=1000))
				}, error = function(e){
				glm(formula, data, family="poisson")
				}, warning=function(w) {				
				glm(formula, data, family="poisson")
				})

				df = coef(summary(fit))[coef,,drop=FALSE]
		}
	}

	df
}








#' Bar plot of cell compositions
#'
#' Bar plot of cell compositions
#'
#' @param obj \code{SingleCellExperiment} from \code{aggregateToPseudoBulk}
#' @param col array of colors.  If missing, use default colors.  If \code{names(col)} is the same as \code{arrayNames(obj)}, then colors will be assigned by assay name#' 
#' @param width specify width of bars
#'
#' @importFrom variancePartition plotPercentBars ggColorHue
#' @export
plotCellComposition = function(obj, col, width=NULL){

  # extract cell counts and other meta-data
  df_cellCount = do.call(rbind, int_colData(obj)$n_cells)

  df_frac = apply(df_cellCount, 1, function(x) x / sum(x))  
  df = as.data.frame(t(df_frac))

  if( missing(col) ){
  	col = ggColorHue(ncol(df))
  }else if( identical(sort(names(col)), sort(colnames(df))) ){
  	col = col[colnames(df)]  	
  }else if( length(col) < ncol(df) ){
  	stop("Too few colors specified: ", length(col), ' < ', ncol(df) )
  }

  plotPercentBars( df, col = col, width=width ) + ylab("Cell percentage")
}




















