

#' Compute logit ratios and their precisions
#' 
#' Compute logit ratios and their precisions given a matrix of counts
#' 
#' @param data matrix of counts where rows are samples and columns are categories, or object that this information can be extracted from
#' @param pc pseudocount
#' 
#' @return \code{EList} defined by \link{limma} storing logit ratio in \code{E} and precisions in \code{weights}
#' 
#' @details 
#' Let \eqn{c1 \sim \Gamma(\alpha, 1)}, \eqn{c2 \sim \Gamma(\beta, 1)} and \eqn{x = \frac{c1}{c2}}, then \eqn{x \sim B(\alpha, \beta)}, the expected value of \eqn{x} is \eqn{E[x] = \frac{\alpha}{\beta}}.  The expected log is \eqn{E[\log(x)] = \psi(\alpha) - \psi(\alpha - \beta)}, where \eqn{\psi()} is the \code{digamma()} function.  
#' 
#' The variance of the logit transformed fraction is \eqn{var[log(\frac{x}{1-x})] = \psi_1(\alpha) + \psi_1(\beta)}, where \eqn{\psi_1()} is the \code{trigamma()} function. For sufficiently large values of \eqn{\alpha} and \eqn{\beta}, this is well approximated by \eqn{1/\alpha + 1/\beta}.
#' 
#' The counts producing these ratios are discrete.  But approximating the counts by \eqn{\Gamma(\lambda, 1)} instead of \eqn{Pois(\lambda)} matches the first two moments exactly.
#' 
#' See \url{https://en.wikipedia.org/wiki/Beta_distribution#Moments_of_logarithmically_transformed_random_variables}
#' 
#' Since \code{log()} and \code{trigamma()} are not defined at zero, add a pseudocount of 0.25.
#' @export
setGeneric("logitRatios", 
	function( data, pc = 0.25){

	standardGeneric("logitRatios")
})


#' @export
#' @rdname logitRatios
#' @aliases logitRatios,matrix-method
setMethod("logitRatios", "matrix",
	function( data, pc = 0.25){

	# add pseudocount
	countMatrix = data + pc

	# compute logit of count ratio
	# countMatrix[1,1] / sum(countMatrix[1,])
	# boot::inv.logit(logitRatios[1,1])
	logitRatios = apply(countMatrix, 1, function(counts) log(counts) - log(sum(counts) - counts))

	# compute variance of logitRatios based on counts
	# https://stats.stackexchange.com/questions/118403/logit-standard-error
	# logitVar = apply(countMatrix, 1, function(counts) 1/counts + 1/(sum(counts) - counts))
	logitVar = apply(countMatrix, 1, function(counts) trigamma(counts) + trigamma(sum(counts) - counts))

	# convert to precision weights
	weights = 1 / logitVar
	
	new("EList", list(E = logitRatios, weights = weights))
})


#' @export
#' @rdname logitRatios
#' @aliases logitRatios,SingleCellExperiment-method
setMethod("logitRatios", "SingleCellExperiment",
	function( data, pc = 0.25){

	# extract cell counts from SingleCellExperiment
	countMatrix = do.call(rbind, int_colData(data)$n_cells)

	logitRatios(countMatrix, pc)
})







# library(boot)

# c1 = 10
# c2 = 7

# C1 = rpois(1000000, c1) + .25
# C2 = rpois(1000000, c2)+ .25
# x = rnorm(100000)

# ratio = C1 / (C1+C2)

# mean(logit(ratio))
# var(logit(ratio))

# # https://stats.stackexchange.com/questions/118403/logit-standard-error
# 1/mean(C1) + 1/mean(C2)

# trigamma(c1) + trigamma(c2)




# x = rpois(100000, 10)
# mean(x)
# var(x)


# x = rgamma(100000, 10, 1)
# mean(x)
# var(x)


# x = seq(1, 30)
# plot(x, dpois(round(x), 10))

# lines(x, dgamma(x, 10, 1), col="red")
















