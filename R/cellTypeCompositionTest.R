

#' Test difference in cell type composition
#'
#' Test difference in cell type composition using a regression models
#'
#' @param obj \code{SingleCellExperiment} from \code{aggregateToPseudoBulk}
#' @param formula formula indicating variables to include in the regression
#' @param coef which coefficient to test
#' @param method which regression model to use
#'
#' @details Evaluates regression model where cell fraction is response using one of the 6 specified methods.  Methods can be fit with random effects, except for 'betabinomial'
#'
#' @importFrom variancePartition dream voomWithDreamWeights eBayes topTable
#' @importFrom edgeR DGEList
#' @export
cellTypeCompositionTest = function( obj, formula, coef, method = c( "nb", "binomial", "betabinomial", "lm", "lmlog", "poisson", "dream")  ){

	stopifnot(is( obj, "SingleCellExperiment") )
	method = match.arg(method)

	# extract cell counts from SingleCellExperiment
	countMatrix = do.call(rbind, int_colData(obj)$n_cells)

	# merge counts with metadata to ensure row ordering
	data = merge(countMatrix, colData(obj), by="row.names")
	rownames(data) = data$Row.names
	data = data[,-1]

	# extract columns now that order is ensured
	countMatrix = data[,seq_len(ncol(countMatrix)),drop=FALSE]
	data = data[,-seq_len(ncol(countMatrix)),drop=FALSE]

	if( method == "dream"){

		dge = DGEList(t(countMatrix))
		vobj = voomWithDreamWeights(dge, ~1, data)
		fit = dream(vobj, formula, data)
		fit = eBayes(fit)

		tab = topTable(fit, coef=coef, number=Inf, sort.by="none")
		tab = with(tab, data.frame(assay = rownames(tab), Estimate=logFC, se=logFC/t, zstat=t, pValue=P.Value, p.adj=adj.P.Val))

		}else{
			tab = testComposition(countMatrix, formula, data, coef, eval="test", method)
		}

		tab
}

#' Variance partitioning of cell type composition
#'
#' Variance partitioning of cell type composition using a regression models
#'
#' @param obj \code{SingleCellExperiment} from \code{aggregateToPseudoBulk}
#' @param formula formula indicating variables to include in the regression
#' @param method which regression model to use
#'
#' @details Evaluates regression model where cell fraction is response using one of the 6 specified methods.  Methods can be fit with random effects, except for 'betabinomial'
#'
#' @export
cellTypeCompositionVarPart = function( obj, formula, method = c( "nb", "binomial", "betabinomial", "lm", "lmlog", "poisson")  ){

	stopifnot(is( obj, "SingleCellExperiment") )
	method = match.arg(method)

	# extract cell counts from SingleCellExperiment
	countMatrix = do.call(rbind, int_colData(obj)$n_cells)

	# merge counts with metadata to ensure row ordering
	data = merge(countMatrix, colData(obj), by="row.names")
	rownames(data) = data$Row.names
	data = data[,-1]

	# extract columns now that order is ensured
	countMatrix = data[,seq_len(ncol(countMatrix)),drop=FALSE]
	data = data[,-seq_len(ncol(countMatrix)),drop=FALSE]

	testComposition(countMatrix, formula, data, coef=NULL, eval="vp", method)
}




# @param eval evaluate either hypothesis test 'test' or variance fractions 'vp'
testComposition = function( countMatrix, formula, data, coef, eval = c("test", "vp"), method = c( "nb", "binomial", "betabinomial", "lm", "lmlog", "poisson") ){

	method = match.arg(method)
	eval = match.arg(eval)

	# extract cell counts and other meta-data
	data = data.frame(data, countMatrix, TotalCells = rowSums(countMatrix), check.names=FALSE)

  	# loop thru assays and perform test
	res = lapply( colnames(countMatrix), function(k){

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
		res = regModel(form, data, coef, eval, method)

		data.frame(assay = k, res)
	})

	# combine results
	df = as.data.frame(do.call(rbind, res))

	if( eval == "test"){
		colnames(df)[3:5] = c("se", "zstat", "pValue")
		df$p.adj = p.adjust(df$pValue, "fdr")
	}

	df
}

#' @importFrom lme4 findbars
.isMixedModelFormula = function(formula){	
    !is.null(findbars(as.formula(formula)))
}

#' @importFrom aod betabin
#' @importFrom MASS glm.nb
#' @importFrom lme4 glmer glmer.nb lmerControl glmerControl .makeCC
#' @importFrom lmerTest lmer
#' @importFrom stats lm glm glm.control
regModel = function(formula, data, coef, eval = c("test", "vp"), method = c( "nb", "binomial", "betabinomial", "lm", "lmlog", "poisson") ){

	method = match.arg(method)
	eval = match.arg(eval)

	formula = as.formula(formula)

	if( .isMixedModelFormula(formula) ){

		# mixed model
		if( method %in% c("binomial",  "poisson") ){
			control = glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
			fit = glmer(formula, data, family=method, control=control)
			df_coef = coef(summary(fit))

		}else if(method == "betabinomial"){
			stop("'betabinomial' model can't accept random effects")

		}else if(method %in% c("lm", "lmlog") ){
			control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
	
			fit = lmer(formula, data, control = control)			
			colids = c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
			df_coef = coef(summary(fit))[,colids,drop=FALSE]

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
			df_coef = coef(summary(fit))
		}
	}else{
		# only fixed effects
		if( method %in% c("binomial",  "poisson") ){
			fit = glm(formula, data, family=method)
			df_coef = coef(summary(fit))

		}else if(method == "betabinomial"){
			fit = betabin(formula, ~ 1, data)
			df_coef = summary(fit)@Coef[coef,,drop=FALSE]

		}else if(method %in% c("lm", "lmlog") ){
			fit = lm(formula, data)
			df_coef = coef(summary(fit))

		}else if(method == "nb" ){
			# if NB fit fails, use poisson model
			fit = tryCatch({
				glm.nb(formula, data, control=glm.control(maxit=1000))
				}, error = function(e){
				glm(formula, data, family="poisson")
				}, warning=function(w) {				
				glm(formula, data, family="poisson")
				})

				df_coef = coef(summary(fit))
		}
	}

	# which result to return
	if( eval == "test"){

		if( ! coef %in% rownames(df_coef) ){
			stop(paste0("coef is not found in model fit....\n coef: ", coef, "\n Coefficients: ",paste(rownames(df_coef), collapse=', '), "\n Note: only fixed effects can be tested" ))
		}

		res = df_coef[coef,,drop=FALSE]
	}else{
		res = calcVarPart( fit )
		res = as.data.frame(t(res))
	}

	res
}



# library(edgeR)
# dge = DGEList(df_cellCount)
# vobj = voomWithDreamWeights(dge, ~1, data.frame(mu = rep(1,8)), plot=TRUE)
















