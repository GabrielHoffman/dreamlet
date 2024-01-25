#' Multivariate outlier detection
#'
#' Detect multivariante outliers using Mahalanobis distance using mean and covariance estimated either with standard or robust methods.
#'
#' @param data matrix of data
#' @param robust use robust methods, defaults to \code{FALSE}
#' @param ... arguments passed to \code{MASS::cov.rob()}
#'
#' @details
#' The distance follow a chisq distrubtion under the null with standard method for mean and covariance.  It is approximate if the robust method is used.  So use \code{qchisq(p = 0.999 , df = k)} to get cutoff to keep 99.9\% of samples under the null for data with \code{k=2} columns.
#'
#' @return \code{data.frame} storing chisq and z-score for each entry indicating deviation from the mean.  The z-score is computed by evaluating the p-value of chisq statistic and converting it into a z-score
#'
#' @examples
#' data <- matrix(rnorm(200), 100, 2)
#'
#' res <- outlier(data)
#'
#' res[1:4,]
#' @importFrom MASS cov.rob
#' @importFrom stats cov mahalanobis pchisq prcomp qnorm
#' @export
outlier <- function(data, robust = FALSE, ...) {

  if (robust) {
    res <- cov.rob(data, ...)

    mu <- res$center
    C <- res$cov
  } else {
    # means
    mu <- colMeans2(data, useNames=FALSE)
    # covariance
    C <- cov(data)
  }

  # Mahalanobis distance
  d = mahalanobis(data, mu, C)

  # P-values based on d ~ chisq(k)
  pOut <- pchisq(d, ncol(data), lower.tail=FALSE)
  z <- qnorm(pOut/2, lower.tail=FALSE)

  df = data.frame(chisq = d, z = z)
  rownames(df) = rownames(data)
  df
}


 