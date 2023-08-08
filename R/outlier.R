#' Multivariate outlier detection
#'
#' Detect multivariante outliers using Mahalanobis distance using mean and covariance estimated either with standard or robust methods.
#'
#' @param data matrix of data
#' @param robust use robust methods
#' @param ... arguments passed to \code{MASS::cov.rob()}
#'
#' @details
#' The distance follow a chisq distrubtion under the null with standard method for mean and covariance.  It is approximate of robust method is used.  So use \code{qchisq(p = 0.999 , df = 2)} to get cutoff to keep 99.9\% of samples under the null for data with 2 columns.
#'
#' @return z-score for each entry indicating deviation from the mean
#'
#' @examples
#' data <- matrix(rnorm(200), 100, 2)
#'
#' res <- outlier(data)
#'
#' res[1:4]
#' @importFrom MASS cov.rob
#' @importFrom stats cov mahalanobis
#' @export
outlier <- function(data, robust = TRUE, ...) {
  if (robust) {
    res <- cov.rob(data, ...)

    mu <- res$center
    C <- res$cov
  } else {
    # means
    mu <- colMeans(data)
    # covariance
    C <- cov(data)
  }

  # Mahalanobis distance
  mahalanobis(data, mu, C)
}
