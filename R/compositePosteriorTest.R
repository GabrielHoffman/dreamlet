# Gabriel Hoffman
# Dec 4, 2023
#
# Perform composite test on results from mashr


#' Perform composite test on results from mashr
#'
#' Perform composite test evaluating the specificity of an effect.  Evalute the posterior probability that an a non-zero effect present in _all_ or _at least one_ condition in the inclusion set, but _no conditions_ in the exclusion set.
#'
#' @param x \code{"dreamlet_mash_result"} from \code{run_mash()}
#' @param include array of conditions in the inclusion set
#' @param exclude array of conditions in the exclusion set. Defaults to \code{NULL} for no exclusion
#' @param test evaluate the posterior probability of a non-zero effect in \code{"at least 1"} or \code{"all"} conditions
#'
#' @description The posterior probabilities for all genes and conditions is obtained as \code{1-lFSR}.  Let \code{prob} be an array storing results for one gene.  The probability that _no_ conditions in the exclusion set are non-zero is \code{prod(1 - prob[exclude])}. The probability that _all_ conditions in the inclusion set are non-zero is \code{prod(prob[include])}. The probability that _at least one_ condition in the inclusion set is non-zero is \code{1 - prod(1 - prob[include])}.  The composite test is the product of the probabilties computed from the inclusion and exclusion sets.
#'
#' @seealso \code{run_mash()}
#' @examples
#' library(muscat)
#' library(mashr)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce[1:100, ],
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#'
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~group_id)
#'
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data
#' res.dl <- dreamlet(res.proc, ~group_id)
#'
#' # run MASH model
#' # This can take 10s of minutes on real data
#' # This small datasets should take ~30s
#' res_mash <- run_mash(res.dl, "group_idstim")
#'
#' # Composite test based on posterior probabilities
#' # to identify effect present in *at least 1* monocyte type
#' # and *NO* T-cell type.
#' include <- c("CD14+ Monocytes", "FCGR3A+ Monocytes")
#' exclude <- c("CD4 T cells", "CD8 T cells")
#'
#' # Perform composite test
#' prob <- compositePosteriorTest(res_mash, include, exclude)
#'
#' # examine the lFSR for top gene
#' get_lfsr(res_mash$model)[which.max(prob), , drop = FALSE]
#' 
#' # Test if *all* cell types have non-zero effect
#' prob <- compositePosteriorTest(res_mash, assayNames(res.dl))
#
#' @export
compositePosteriorTest <- function(x, include, exclude = NULL, test = c("at least 1", "all")) {
  test <- match.arg(test)

  if( is(x, "dreamlet_mash_result") ){
    # get probability from lFSR
    prob <- 1 - get_lfsr(x$model)
  }else if(is(x, "data.frame") | is(x, "matrix")){
    prob = x
  }else{
    stop("data type of x not recognized")
  }

  .compositePosteriorTest(prob, include, exclude, test)
}

# Given matrix of posterior probabilities with genes on rows and columns as conditions, compute composite probability from include vs exclude set.
.compositePosteriorTest <- function(prob, include, exclude=NULL, test = c("at least 1", "all")) {
  test <- match.arg(test)

  # if probability is NA, set it to zero
  prob[is.na(prob)] <- 0

  stopifnot(all(include %in% colnames(prob)))
  stopifnot(all(exclude %in% colnames(prob)))

  # probability that *NO* cell types have a non-zero effect
  prob_excl <- apply(1 - prob[, exclude, drop = FALSE], 1, prod, na.rm = TRUE)

  if (test == "at least 1") {
    # Probability at least 1, (i.e. probability that not none)
    prob_incl <- 1 - apply(1 - prob[, include, drop = FALSE], 1, prod, na.rm = TRUE)
  } else if (test == "all") {
    # for each gene
    # probability that *ALL* cell types have a non-zero effect
    prob_incl <- apply(prob[, include, drop = FALSE], 1, prod, na.rm = TRUE)
  }

  prob_incl * prob_excl
}




