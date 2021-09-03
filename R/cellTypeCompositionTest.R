

#' Test difference in cell type composition
#'
#' Test difference in cell type composition using a regression models
#'
#' @param obj \code{SingleCellExperiment} from \code{aggregateToPseudoBulk}
#' @param formula formula indicating variables to include in the regression
#' @param coef which coefficient to test
#' @param method which regression model to use
#'
#' @importFrom aod betabin
#' @importFrom MASS glm.nb
#'
#' @export
cellTypeCompositionTest = function( obj, formula, coef, method = c("binomial", "betabinomial", "lm", "lmlog", "poisson", "nb") ){

	method = match.arg(method)

	# extract cell counts and other meta-data
	df_cellCount = do.call(rbind, int_colData(obj)$n_cells)
	df_cellCount = cbind(df_cellCount, TotalCells = rowSums(df_cellCount))
	df_cellCount = merge(df_cellCount, colData(pb), by="row.names")

  	# loop thru assays and perform test
	df_fit = lapply( assayNames(pb), function(k){

		cat(k)
		# create formula
		form = paste0("cbind(`", k, "`, TotalCells - `", k, "`) ~ ", as.character(formula)[2])

	    # evaluate regression
	   	res = switch( method,
			"binomial" = {
				fit = glm(as.formula(form), data=df_cellCount, family="binomial")
				coef(summary(fit))[coef,,drop=FALSE]},
			"betabinomial" = { 
				fit = betabin(as.formula(form), ~ 1, data=df_cellCount)
				summary(fit)@Coef[coef,,drop=FALSE]},
			"lm" = {
				form = paste0("`", k, "`  / TotalCells ~ ", as.character(formula)[2])
				fit = lm(as.formula(form), data=df_cellCount)
				coef(summary(fit))[coef,,drop=FALSE]},
			"lmlog" = {
				form = paste0("log(`", k, "`  / TotalCells) ~ ", as.character(formula)[2])
				fit = lm(as.formula(form), data=df_cellCount)
				coef(summary(fit))[coef,,drop=FALSE]},
			"poisson" = {
				form = paste0("`", k, "` ~ offset(log(TotalCells)) + ", as.character(formula)[2])
				fit = glm(as.formula(form), data=df_cellCount, family="poisson")
				coef(summary(fit))[coef,,drop=FALSE]},
			"nb" = {
				form = paste0("`", k, "` ~ offset(log(TotalCells)) + ", as.character(formula)[2])
				# if error or warning in NB, use Poisson
				tryCatch({
					fit = glm.nb(as.formula(form), data=df_cellCount, control=glm.control(maxit=1000))
					coef(summary(fit))[coef,,drop=FALSE]}, 
					error = function(e){
						fit = glm(as.formula(form), data=df_cellCount, family="poisson")
						coef(summary(fit))[coef,,drop=FALSE]
					},
					warning=function(w) {
						fit = glm(as.formula(form), data=df_cellCount, family="poisson")
						coef(summary(fit))[coef,,drop=FALSE]
					})
				}
	      )
		data.frame(assay =k, res)
	})
	df_fit = as.data.frame(do.call(rbind, df_fit))
	colnames(df_fit)[3:5] = c("se", "zstat", "pValue")
	df_fit$p.adj = p.adjust(df_fit$pValue, "fdr")

	df_fit
}
