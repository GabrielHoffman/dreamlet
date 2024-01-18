# Gabriel Hoffman
# Nov 15, 2021

#' Class dreamlet_mash_result
#'
#' Class \code{dreamlet_mash_result}
#'
#' @name dreamlet_mash_result-class
#' @rdname dreamlet_mash_result-class
#' @exportClass dreamlet_mash_result
#' @return \code{dreamlet_mash_result} class
setClass("dreamlet_mash_result", contains = "list")



#' Convert results table to matrix
#'
#' Convert results table to matrix
#'
#' @param tab results table from \code{topTable()}
#' @param col which column to extract
#' @param rn column id storing rownames
#' @param cn column id storing colnames
#'
#' @return matrix storing values of column \code{col} in rows defind by \code{rn} and columns defined by \code{cn}
#'
#' @importFrom Matrix sparseMatrix
tabToMatrix <- function(tab, col, rn = "ID", cn = "assay") {
  # check column names
  if (!rn %in% colnames(tab)) stop("Column not found: ", rn)
  if (!cn %in% colnames(tab)) stop("Column not found: ", cn)

  # convert query column names to factor
  if (!is.factor(tab[[rn]])) tab[[rn]] <- factor(tab[[rn]])
  if (!is.factor(tab[[cn]])) tab[[cn]] <- factor(tab[[cn]])

  # extract row and column names for resulting matrix
  rnlvl <- levels(tab[[rn]])
  cnlvl <- levels(tab[[cn]])

  # get indeces
  i <- match(tab[[rn]], rnlvl)
  j <- match(tab[[cn]], cnlvl)

  # convert row,col,value to sparse matrix
  # empty entries are set to 0
  M <- sparseMatrix(i, j,
    x = tab[[col]],
    dims = c(length(rnlvl), length(cnlvl)),
    dimnames = list(rnlvl, cnlvl)
  )

  # convert to real matrix
  data <- as.matrix(M)

  # replace 0's with NA
  # since
  data[data == 0] <- NA

  data
}


#' Run mash analysis on dreamlet results
#'
#' Run mash analysis on dreamlet results
#'
#' @param fit result from \code{dreamlet()}
#' @param coefList coefficient to be analyzed
#'
#' @details
#' Apply \href{https://cran.r-project.org/web/packages/mashr/index.html}{mashr} analysis \insertCite{urbut2019flexible}{dreamlet} on the joint set of coefficients for each gene and cell type.  \code{mashr} is a Bayesian statistical method that borrows strength across tests (i.e. genes and cell types) by learning the distribution of non-zero effects based the obesrved logFC and standard errors.  The method then estimates the posterior distributions of each coefficient based on the observed value and the genome-wide emprical distribution.
#'
#' \code{mashr} has been previously applied to differential expression in \href{https://www.gtexportal.org}{GTEx} data using multiple tissues from the same set of donors \insertCite{oliva2020impact}{dreamlet}.
#'
#' In single cell data, a given gene is often not sufficiently expressed in all cell types.  So it is not evaluated in a subsets of cell types, and its coefficient value is \code{NA}. Since mashr assumes coefficients and standard errors for every gene and cell type pair, entries with these missing values are set to have \code{coef = 0}, and \code{se = 1e6}.  The output of mashr is then modified to set the corresponding values to \code{NA}, to avoid nonsensical results downstream.
#'
#' @return a list storing the \code{mashr} model as \code{model} and the original coefficients as \code{logFC.original}
#'
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
#' # extract statistics from mashr model
#' # NA values indicate genes not sufficiently expressed
#' # in a given cell type
#'
#' # original logFC
#' head(res_mash$logFC.original)
#'
#' # posterior mean for logFC
#' head(get_pm(res_mash$model))
#'
#' # how many gene-by-celltype tests are significant
#' # i.e.  if a gene is significant in 2 celltypes, it is counted twice
#' table(get_lfsr(res_mash$model) < 0.05, useNA = "ifany")
#'
#' # how many genes are significant in at least one cell type
#' table(apply(get_lfsr(res_mash$model), 1, min, na.rm = TRUE) < 0.05)
#'
#' # how many genes are significant in each cell type
#' apply(get_lfsr(res_mash$model), 2, function(x) sum(x < 0.05, na.rm = TRUE))
#'
#' # examine top set of genes
#' # which genes are significant in at least 1 cell type
#' sort(names(get_significant_results(res_mash$model)))[1:10]
#'
#' # Lets examine ENO1
#' # There is a lot of variation in the raw logFC
#' res_mash$logFC.original["ENO1", ]
#'
#' # posterior mean after borrowing across cell type and genes
#' get_pm(res_mash$model)["ENO1", ]
#'
#' # forest plot based on mashr results
#' plotForest(res_mash, "ENO1")
#'
#' # volcano plot based on mashr results
#' # yaxis uses local false sign rate (lfsr)
#' plotVolcano(res_mash)
#'
#' # Comment out to reduce package runtime
#' # gene set analysis using mashr results
#' # library(zenith)
#' # go.gs = get_GeneOntology("CC", to="SYMBOL")
#' # df_gs = zenith_gsa(res_mash, go.gs)
#'
#' # Heatmap of results
#' # plotZenithResults(df_gs, 2, 1)
#'
#' @references{
#' \insertAllCited{}
#' }
#'
#' @seealso \code{mashr::mash_estimate_corr_em()}, \code{mashr::cov_canonical}, \code{mashr::mash_set_data}
#' @importFrom mashr mash_set_data cov_canonical mash_estimate_corr_em
#' @importFrom MatrixGenerics colVars
#' @export
run_mash <- function(fit, coefList) {
  if (!is(fit, "dreamletResult")) {
    stop("fit must be of class dreamletResult")
  }

  if (!coefList %in% coefNames(fit)) {
    stop("coef not found in coefNames(fit): ", coefList)
  }

  tab <- lapply(coefList, function(coef) {
    # get results for each gene and cell type
    tab <- topTable(fit, coef = coef, Inf)
    tab$coef <- coef
    tab
  })
  tab <- do.call(rbind, tab)

  if (length(coefList) > 1) {
    tab$assay <- paste(tab$assay, tab$coef, sep = ".")
    lvls <- c(vapply(assayNames(fit), function(x) paste(x, coefList, sep = "."), character(0)))
  } else {
    lvls <- assayNames(fit)
  }

  # compute standard error from t-stat and logFC
  tab$se <- tab$logFC / tab$t

  # sort tab results based on assayNames(fit)
  tab$assay <- factor(tab$assay, lvls)

  # convert to matricies
  B <- tabToMatrix(tab, "logFC")
  S <- tabToMatrix(tab, "se")

  # only keep columns with variance in logFC
  cv <- colVars(B, na.rm = TRUE, useNames = TRUE)
  keep <- (cv > 0) & !is.na(cv)
  B <- B[, keep, drop = FALSE]
  S <- S[, keep, drop = FALSE]

  # run mashr on these matricies
  #-----------------------------

  # set up
  # NA's are replaced with beta = 0 with se = 1e6
  data <- mash_set_data(B, S)

  # estimate model parameters
  U.c <- cov_canonical(data)

  # Estimate correlation structure
  V.em <- mash_estimate_corr_em(data, U.c, details = TRUE)

  # copy model to drop NA terms
  model <- V.em$mash.model

  # B has the same ordering as these, so replace corresponding elements with NA
  # this revents non-sensiccal results for coefficients that were originally NA
  idx <- which(is.na(B))
  model$result$PosteriorMean[idx] <- NA
  model$result$PosteriorSD[idx] <- NA
  model$result$NegativeProb[idx] <- NA
  model$result$lfsr[idx] <- NA

  # format results as new object
  new("dreamlet_mash_result", list(model = model, logFC.original = B, coefList = coefList))
}
