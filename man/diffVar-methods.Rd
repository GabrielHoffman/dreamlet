% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffVar.R
\name{diffVar,dreamletResult-method}
\alias{diffVar,dreamletResult-method}
\alias{diffVar,dreamletResult,dreamletResult-method}
\title{Test differential variance}
\usage{
\S4method{diffVar}{dreamletResult}(
  fit,
  method = c("AD", "SQ"),
  scale = c("leverage", "none"),
  BPPARAM = SerialParam(),
  ...
)
}
\arguments{
\item{fit}{model fit from \code{dream()}}

\item{method}{transform the residuals using absolute deviation ("AD") or squared deviation ("SQ").}

\item{scale}{scale each observation by "leverage", or no scaling ("none")}

\item{BPPARAM}{parameters for parallel evaluation}

\item{...}{other parameters passed to \code{dream()}}
}
\description{
Test the association between a covariate of interest and the response's deviation from expectation.
}
\details{
This method performs a test of differential variance between two subsets of the data, in a way that generalizes to multiple categories, continuous variables and metrics of spread beyond variance.  For the two category test, this method is simular to Levene's test.  This model was adapted from Phipson, et al (2014), extended to linear mixed models, and adapted to be compatible with \code{variancePartition::dream()} and \code{dreamlet::dreamlet()}.

This method is composed of multiple steps where 1) a typical linear (mixed) model is fit with \code{dreamlet()}, 2) residuals are computed and transformed based on an absolute value or squaring transform, 3) a second regression is performed with \code{dreamlet()} to test if a variable is associated with increased deviation from expectation.  Both regression take advantage of the \code{dreamlet()} linear (mixed) modelling framework followed by empirical Bayes shrinkage that extends the \code{limma::voom()} framework.

Note that \code{diffVar()} takes the results of the first regression as a parameter to use as a starting point.
}
\examples{
library(muscat)
library(SingleCellExperiment)

data(example_sce)

# create pseudobulk for each sample and cell cluster
pb <- aggregateToPseudoBulk(example_sce,
  assay = "counts",
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  verbose = FALSE
)

# voom-style normalization
res.proc <- processAssays(pb, ~group_id)

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)

# Differential variance analysis
# result is a dreamlet fit
res.dvar <- diffVar(res.dl)

# Examine results
res.dvar

# Examine details for each assay
details(res.dvar)

# show coefficients estimated for each cell type
coefNames(res.dvar)

# extract results using limma-style syntax
# combines all cell types together
# adj.P.Val gives study-wide FDR
topTable(res.dvar, coef = "group_idstim", number = 3)

# Plot top hit to see differential variance
# Note that this is a toy example with only 4 samples
cellType <- "CD4 T cells"
gene <- "DYNLRB1"

y <- res.proc[[cellType]]$E[gene, ]
x <- colData(res.proc)$group_id

boxplot(y ~ x,
  xlab = "Stimulation status",
  ylab = "Gene expression",
  main = paste(cellType, gene)
)
#
}
\references{
{
  \insertRef{phipson2014diffvar}{variancePartition}
}
}
\seealso{
\code{variancePartition::diffVar()}

\code{variancePartition::diffVar()}, \code{missMethyl::diffVar()}
}
