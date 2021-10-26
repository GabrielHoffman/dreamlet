# Gabriel Hoffman
# Oct 25, 2021

# create class to store resilt of compositionTest
# @importFrom variancePartition MArrayLM2
# setClass("MArrayLM2_ilr", representation("MArrayLM2", V.back = "matrix"))

#' Class MArrayLM2_ilr
#'
#' Class \code{MArrayLM2_ilr} 
#'
#' @name MArrayLM2_ilr-class
#' @rdname MArrayLM2_ilr-class
#' @exportClass MArrayLM2_ilr
setClass("MArrayLM2_ilr", contains="MArrayLM" )

setIs("MArrayLM2_ilr","LargeDataObject")

setAs(from='MArrayLM2', to='MArrayLM2_ilr', function(from){
	structure(from, class="MArrayLM2_ilr")
	})

setAs(from='MArrayLM', to='MArrayLM2_ilr', function(from){
	structure(from, class="MArrayLM2_ilr")
	})



# how to retain original names after back transform???


#' Regression testing for composition changes 
#' 
#' Regression testing for composition changes using linear (mixed) models using log ratio of counts
#' 
#' @param counts count data
#' @param formula regression formula for independent variables
#' @param info data.frame storing variables in formula
#' @param pseudocount added to counts to avoid issues with zeros
#' @param useWeights default: TRUE. If TRUE use precision weights, else ignore weights
#' @param BPPARAM parameters for parallel evaluation
#' 
#' @importFrom compositions ilrBase
#' @importFrom variancePartition dream eBayes
#' @importFrom methods is as
#' @export
compositionTest = function(counts, formula, info, pseudocount=.1, useWeights=TRUE, BPPARAM=SerialParam()){

	# transform to rotated log ratios of counts
	# cobj = t(as.data.frame(ilr(counts + pseudocount)))
	# cobj = new("EList", list(E = cobj, weights = cobj))
	# cobj$weights[] =1

	# eval ILR transform with precision weights
	cobj = ilrWithPrecisionWeights( counts, pseudocount=pseudocount )

	if( ! useWeights ){
		cobj$weights[] = 1
	}

	# get matrix to perform back-transfrom from ilr -> clr
	V.back = ilrBase(z=t(cobj$E))
	colnames(V.back) = paste0("ilr_", seq_len(ncol(V.back)))
	rownames(V.back) = colnames(counts)

	# fit regression models
	fit = suppressMessages(dream( cobj, formula, info, BPPARAM=BPPARAM))
	fit = eBayes(fit)

	# return object with model fit and back-transform matrix
	res = as(fit, "MArrayLM2_ilr")

	# store V.back to backtransform later
	res$V.back = V.back

	# store type information about current fit
	res$original = is(fit)

	res
}

# need to return variances, rathern than variance fractions
# @export
# compositionVP = function(counts, formula, info, pseudocount=.1, BPPARAM=SerialParam()){

# 	# eval ILR transform with precision weights
# 	cobj = ilrWithPrecisionWeights( counts, pseudocount=pseudocount )

# 	# fit regression models
# 	vc.ilr = fitVarPartModel( cobj, formula, info, BPPARAM=BPPARAM, fxn=function(fit){
# 		calcVarPart(fit, returnFraction=FALSE)
# 	})
# 	vc.ilr = do.call(rbind, vc.ilr)

# 	# get matrix to perform back-transfrom from ilr -> clr
# 	V.back = ilrBase(z=t(vc.ilr))

# 	# Does ilr -> clr transform for variances work for variance fractions?
# 	vp.backtransform = diag(V.back %*% tcrossprod(vc.ilr, V.back))

# 	# enforce closure so components sum to 1
# 	apply(vp.backtransform, 1, function(x) x / sum(x))
# }




#' Apply isometric log ratio (ilr) transform and compute precision weights
#'
#' Apply isometric log ratio (ilr) transform and compute precision weights
#'
#' @param counts data matrix of counts
#' @param pseudocount pseudocounts to avoid problems with zero counts
#'
#' @description Compute isometric log ratio (ilr) transform of counts, and their observation-level variances using the asymptotic normal approximation of Egozcue, et al (2020).  Convert observation-level variances to precisions in \code{Elist} for downstream use in \code{dream()}.
#'
#' @references{
#'   \insertRef{egozcue2020some}{dreamlet}
#' }
#'
#' @importFrom compositions ilrBase
#' @importFrom methods new
#' @import Rdpack
#' @export
ilrWithPrecisionWeights = function(counts, pseudocount = .1){

	# get ILR transformation matrix
	V = ilrBase(D=ncol(counts))

	# Apply ILR transformation
	# y = t(ilr(counts + pc))
	y = t(log( counts + pseudocount ) %*% V	)
	rownames(y) = paste0("ilr_", seq_len(nrow(y)))

	# compute observation-level theoretical variances
	y_vars = apply(counts + pseudocount, 1, function(x){
		frac = x / sum(x)
		diag(crossprod(V, diag(1/frac)) %*% V) / sum(x)
		})

	if( !is.matrix(y_vars) ){
		y_vars = t(as.matrix(y_vars))
	}
	
	# combine into EList object
	new("EList", list(E = y, weights = 1/y_vars))
}

#' toptable for MArrayLM2_ilr
#'
#' toptable for MArrayLM2_ilr
#'
#' @param fit fit
#' @param coef coef
#' @param number number
#' @param genelist genelist
#' @param adjust.method adjust.method
#' @param sort.by sort.by
#' @param resort.by resort.by
#' @param p.value p.value
#' @param lfc lfc
#' @param confint confint
#'
#' @return results of toptable
#' @export
#' @importFrom stats qnorm
#' @import limma
#' @rdname toptable-method
#' @aliases toptable,MArrayLM2_ilr-method
#' @importFrom compositions ilrBase
#' @importFrom variancePartition topTable
#' @importFrom stats p.adjust pnorm
#' @importFrom methods as
#' @export
setMethod("topTable", "MArrayLM2_ilr", 
	function(fit,
       coef = NULL,
       number = 10,
       genelist = fit$genes,
       adjust.method = "BH",
       sort.by = "B",
       resort.by = NULL,
       p.value = 1,
       lfc = 0,
       confint = FALSE){

	# extract results using standard topTable
	# then back-transform the results from ilr -> clr space	
	if( fit$original[1] == "MArrayLM"){
		fit.orig = as(fit, "MArrayLM")
	}else{
		fit.orig = as(fit, "MArrayLM2")
	}

	tab = topTable(fit.orig, coef=coef, number=number, genelist=genelist, adjust.method=adjust.method, sort.by=sort.by, resort.by=resort.by, p.value=p.value, lfc=lfc, confint=confint)

	# extract coefficients
	beta_irl = tab$logFC

	# back-transform is motivated by 
	# ilrInv(coef(fit)[,coef], orig = counts )

	# create back-transform matrix: ilr -> clr
	# fit$V.back

	# apply back-transform to coefficients
	beta = tcrossprod(beta_irl, fit$V.back)	

	# back-transform variance
	se = with(tab, logFC/t)
	beta.se = sqrt(diag(fit$V.back %*% tcrossprod(se^2, fit$V.back)))

	# compute t-statistic on clr scale
	tstat = c(beta / beta.se)
	p.value = 2*pnorm(abs(tstat), lower.tail=FALSE)

	#  logFC     AveExpr           t   P.Value adj.P.Val       z.std
	data.frame( 
				logFC 		= beta.se,
				t 			= tstat,
				P.Value 	= p.value,
				adj.P.Val 	= p.adjust(p.value, "fdr"))
})
