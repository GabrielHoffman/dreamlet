#' Multivariate outlier detection
#'
#' Detect multivariante outliers using Mahalanobis distance using mean and covariance estimated either with standard or robust methods.
#'
#' @param data matrix of data
#' @param robust use robust covariance method, defaults to \code{FALSE}
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


 
#' Outlier analysis for each assay
#' 
#' Compute outlier score for each sample in each assay using \code{outlier()}
#' 
#' @param object \code{dreamletProcessedData} from \code{processAssays()}
#' @param assays assays / cell types to analyze
#' @param nPC number of PCs to uses for outlier score with \code{outlier()}
#' @param robust use robust covariance method, defaults to \code{FALSE}
#' @param ... arguments passed to \code{MASS::cov.rob()}
#' 
#' @return
#' \itemize{
#'  \item{\code{ID}:}{sample identifier}
#'  \item{\code{assay}:}{specify assay}
#'  \item{\code{PCs}:}{principal components}
#'  \item{\code{chisq}:}{mahalanobis distance that is distributed as chisq(k) k = nPC if the data is multivariate gaussian}
#'  \item{\code{z}:}{z-score corresponding to the chisq distance}
#' }
#'
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce,
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#'
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~group_id)
#' 
#' # Compute PCs and outlier scores
#' outlierByAssay( res.proc, c("B cells", "CD14+ Monocytes"))
#
#' @importFrom irlba prcomp_irlba
#' @seealso \code{outlier()}
#' @export
outlierByAssay = function(object, assays, nPC=2, robust = FALSE, ...){

  stopifnot(is(object, "dreamletProcessedData"))

   df = lapply(assays, function(id){
    # get normalized expression
    Y <- assay(object, id)$E

    # PCA
    # dcmp <- prcomp(scale(t(Y)))
    if( nPC < min(dim(Y))/2 ){
      dcmp = prcomp_irlba(t(Y), n = nPC, center=TRUE, scale.=TRUE)
    }else{
      dcmp = prcomp(t(Y), center=TRUE, scale.=TRUE)
    }

    # outlier analysis on first 2 PCs
    df_outlier <- outlier(dcmp$x[,seq(nPC)] * dcmp$sdev[seq(nPC)], robust = robust, ...)
   
    data.frame(ID = colnames(Y), 
      assay = id, 
      dcmp$x, 
      df_outlier) %>%
      as_tibble 
  })
  bind_rows(df)
}

