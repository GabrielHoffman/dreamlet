% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPCA.R
\name{plotPCA}
\alias{plotPCA}
\alias{plotPCA,list-method}
\title{Plot PCA of gene expression for an assay}
\usage{
\S4method{plotPCA}{list}(
  object,
  assays = names(object),
  nPC = 2,
  robust = FALSE,
  ...,
  maxOutlierZ = 20,
  nrow = 2,
  size = 2,
  fdr.cutoff = 0.05
)
}
\arguments{
\item{object}{\code{dreamletProcessedData} from \code{processAssays()} or a \code{list} from \code{residuals()}}

\item{assays}{assays / cell types to analyze}

\item{nPC}{number of PCs to uses for outlier score with \code{outlier()}}

\item{robust}{use robust covariance method, defaults to \code{FALSE}}

\item{...}{arguments passed to \code{MASS::cov.rob()}}

\item{maxOutlierZ}{cap outlier z-scores at this value for plotting to maintain consistent color scale}

\item{nrow}{number of rows in plot}

\item{size}{size passed to \code{geom_point()}}

\item{fdr.cutoff}{FDR cutoff to determine outlier}
}
\description{
Compute PCA of gene expression for an assay, and plot samples coloring by outlier score
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

# PCA to identify outliers
# from normalized expression
plotPCA( res.proc, c("B cells", "CD14+ Monocytes"))

# Run on regression residuals
#-----------------------------

# Regression analysis
fit = dreamlet(res.proc, ~ group_id)

# Extract regression residuals
residsObj = residuals(fit)

# PCA on residuals
plotPCA( residsObj, c("B cells", "CD14+ Monocytes"))
}
\seealso{
\code{outlierByAssay()}
}
