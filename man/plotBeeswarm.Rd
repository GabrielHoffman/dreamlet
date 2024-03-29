% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotBeeswarm.R
\name{plotBeeswarm}
\alias{plotBeeswarm}
\title{Beeswarm plot of effect sizes for each assay}
\usage{
plotBeeswarm(res.dl, coef, fdr.range = 4, assays = assayNames(res.dl))
}
\arguments{
\item{res.dl}{\code{dreamletResult} object from \code{dreamlet()}}

\item{coef}{coefficient name fed to \code{topTable()}}

\item{fdr.range}{range for coloring FDR}

\item{assays}{which assays to plot}
}
\value{
\code{ggplot2} of logFC by assay
}
\description{
Beeswarm plot of effect sizes for each assay, colored by sign and FDR
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

# Beeswarm plot of effect sizes for each assay,
# colored by sign and FDR
plotBeeswarm(res.dl, "group_idstim")

}
