# Gabriel Hoffman
# October 26, 2021

#' Centered log ratio transform
#' 
#' Compute the centered log ratio (clr) transform of a count matrix.
#' 
#' @param counts count data with samples as rows and variables are columns
#' @param pseudocount added to counts to avoid issues with zeros
#' 
#' @seealso compositions::clr
clr = function(counts, pseudocount = 0.5){

	if( ! is.matrix(counts) ){
		counts = matrix(counts, nrow=1)
	}

	log(counts + pseudocount) - rowMeans(log(counts + pseudocount))
}

#' Count ratio uncertainty modeling base linear regression
#' 
#' Count ratio uncertainty modeling base linear regression (crumblr) returns transformed counts and 
#' 
#' @param counts count data with samples as rows and variables are columns
#' @param pseudocount added to counts to avoid issues with zeros
#' 
#' @return  An \code{EList} object with the following components:
#' \itemize{
#'  \item{E:}{numeric matrix of CLR transformed counts}
#'  \item{weights:}{numeric matrix of observation-level inverse-variance weights}
#' }
#' 
#' @description
#' \deqn{\text{clr}_i({\bf \hat p}) = \log(p_i) - \frac{1}{D}\sum_{j=1}^D \log(p_j)}
#'  with sampling variance
#' \deqn{\text{var}[\text{clr}_i({\bf \hat p})] =  \frac{1}{n} \left[ \frac{1}{p_i} - \frac{2}{ D p_i} + \frac{1}{D^2}\sum_{j=1}^D p_j  \right] }
#' 
#' @export
crumblr = function(counts, pseudocount = 0.5){

	D = ncol(counts)
	
	# var_asymp = (1/p - 2/(p*D) + sum(1/p)/D^2) / n
	var_asymp = apply(counts + pseudocount, 1, function(x){
		p = x / sum(x) # fractions
		n = sum(x) # total counts
		(1/p - 2/(p*D) + sum(1/p)/D^2) / n
		})
	
	Y_clr = t(clr(counts, pseudocount))

	new("EList", list(E = Y_clr, weights = 1/var_asymp))
}


