# Gabriel Hoffman
# Jan 16, 2024

#' Meta-analysis across multiple studies
#'
#' Meta-analysis across multiple studies
#'
#' @param x \code{data.frame} rbind'ing results across genes, cell types and datasets
#' @param method meta-analysis method.  Values are fed into \code{metafor::rma()}, except for \code{'RE2C'} which calls \code{remaCor::RE2C()}.
#' @param group colums in \code{x} to group by.  For results from \code{dreamlet::topTable()}, results are aggregrated by gene and cell type (i.e. \code{'ID'} and \code{'assay'}).  If \code{x} is not from this function, this argument allows the function to group results properly
#' @param control passed to \code{rma(..,control)}
#'
#' @details
#' \itemize{
#'  \item{\code{'FE'}: }{fixed effects meta-analysis}
#'  \item{\code{'REML'}: }{random effects meta-analysis}
#'  \item{\code{'RE2C'}: }{joint testing of fixed and random effects}
#' }
#'
#' @examples
#' library(dreamlet)
#' library(muscat)
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
#' # just 'CD14+ Monocytes' for speed
#' res.proc <- processAssays(pb, ~group_id, assays = "CD14+ Monocytes")
#'
#' # dreamlet
#' res.dl <- dreamlet(res.proc, ~group_id)
#'
#' tab1 <- topTable(res.dl, coef = "group_idstim", number = Inf)
#' tab1$Dataset <- "1"
#'
#' # Results from a second cohort
#' # Here, just a copy of the same results for simplicity
#' tab2 <- tab1
#' tab2$Dataset <- "2"
#'
#' # rbind
#' tab_combined <- rbind(tab1, tab2)
#'
#' # Perform fixed effects meta-analysis
#' res <- meta_analysis(tab_combined, method = "FE")
#'
#' res[1:3, ]
#' #
#' @importFrom broom tidy
#' @importFrom metafor rma
#' @importFrom remaCor RE2C
#' @importFrom dplyr select as_tibble group_by_at group_modify
#' @export
meta_analysis <- function(x, method = "FE", group = c("ID", "assay"), control = list(maxiter = 2000)) {
  RE2Cp.twoStep <- term <- type <- NULL

  # evaluate meta-analysis on one gene / assay
  f <- function(beta, se, method) {
    if (method == "RE2C") {
      if (length(beta) > 1) {
        tab <- RE2C(beta, se) %>%
          select(-RE2Cp.twoStep)
      } else {
        tab <- data.frame(stat1 = NA)
      }
    } else {
      tab <- tidy(rma(yi = beta, sei = se, method = method, control = control)) %>%
        select(-term, -type)
    }
    tab$n.studies <- length(beta)
    tab$method <- method
    tab
  }

  # run on each gene and assay
  x %>%
    as_tibble() %>%
    group_by_at(group) %>%
    group_modify(~ f(.x$logFC, .x$logFC / .x$t, method))
}



## res <- tryCatch(
## 	rma(),
## 	error = function(e){NA})
# #if( is.na(res)) return(data.frame())
