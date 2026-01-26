# Plot 2D projection

Plot 2D projection (i.e. UMAP, tSNE) for millions of cells efficiently

## Usage

``` r
plotProjection(
  sce,
  type,
  annotation,
  pointsize = 0,
  pixels = c(512, 512),
  legend.position = "none",
  text = TRUE,
  order
)
```

## Arguments

- sce:

  `SingleCellExperiment`

- type:

  field in `reducedDims(sce)` to plot

- annotation:

  column in `colData(sce)` to annotate each cell

- pointsize:

  Radius of rasterized point. Use `0` for single pixels(fastest).

- pixels:

  Vector with X and Y resolution of the raster, default `c(512,512)`

- legend.position:

  legend.position: the position of legends ("none", "left", "right",
  "bottom", "top", or two-element numeric vector)

- text:

  show `annotation` as text. Default `TRUE`

- order:

  specify order of levels for `annotation`

## Value

ggplot2 plot of the projection

## Details

Uses
[`scattermore::geom_scattermore()`](https://rdrr.io/pkg/scattermore/man/geom_scattermore.html)
to plot millions of points efficiently

## Examples

``` r
library(muscat)
library(SingleCellExperiment)

data(example_sce)

plotProjection(example_sce, "TSNE", "cluster_id", 1)
```
