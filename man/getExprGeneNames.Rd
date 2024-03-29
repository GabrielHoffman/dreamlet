% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getExprGeneNames.R
\name{getExprGeneNames}
\alias{getExprGeneNames}
\title{Get list of expressed genes for each assay}
\usage{
getExprGeneNames(
  sceObj,
  assays = assayNames(sceObj),
  min.cells = 5,
  min.count = 5,
  min.samples = 4,
  min.prop = 0.4,
  min.total.count = 15,
  normalize.method = "TMM"
)
}
\arguments{
\item{sceObj}{SingleCellExperiment object}

\item{assays}{array of assay names to include in analysis. Defaults to \code{assayNames(sceObj)}}

\item{min.cells}{minimum number of observed cells for a sample to be included in the analysis}

\item{min.count}{minimum number of reads for a gene to be considered expressed in a sample.  Passed to \code{edgeR::filterByExpr}}

\item{min.samples}{minimum number of samples passing cutoffs for cell cluster to be retained}

\item{min.prop}{minimum proportion of retained samples with non-zero counts for a gene to be retained}

\item{min.total.count}{minimum total count required per gene for inclusion}

\item{normalize.method}{normalization method to be used by \code{calcNormFactors}}
}
\description{
Get list of expressed genes for each assay using same filters as \code{processAssays()}.
}
\examples{
library(muscat)

data(example_sce)

# create pseudobulk for each sample and cell cluster
pb <- aggregateToPseudoBulk(example_sce,
  assay = "counts",
  sample_id = "sample_id",
  cluster_id = "cluster_id",
  verbose = FALSE
)

# Gene expressed genes for each cell type
geneList = getExprGeneNames(pb)

# Create precision weights for pseudobulk
# By default, weights are set to cell count,
# which is the default in processAssays()
# even when no weights are specified
weightsList <- pbWeights(example_sce,
  sample_id = "sample_id",
  cluster_id = "cluster_id",
  geneList = geneList
)

# voom-style normalization using initial weights
res.proc <- processAssays(pb, ~group_id, weightsList = weightsList)
}
