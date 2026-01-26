# Testing non-linear effects

## Introduction

Typical analysis using regression models assumes a linear affect of the
covariate on the response. Here we consider testing non-linear effects
in the case of 1) continuous and 2) ordered categorical variables.

We demonstrate this feature on a lightly modified analysis of PBMCs from
8 individuals stimulated with interferon-β ([Kang, et al, 2018, Nature
Biotech](https://www.nature.com/articles/nbt.4042)).

## Standard processing

Here is the code from the main vignette:

``` r

library(dreamlet)
library(muscat)
library(ExperimentHub)
library(scater)

# Download data, specifying EH2259 for the Kang, et al study
eh <- ExperimentHub()
sce <- eh[["EH2259"]]

# only keep singlet cells with sufficient reads
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
sce <- sce[, colData(sce)$multiplets == "singlet"]

# compute QC metrics
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

# set variable indicating stimulated (stim) or control (ctrl)
sce$StimStatus <- sce$stim

sce$id <- paste0(sce$StimStatus, sce$ind)

# Create pseudobulk
pb <- aggregateToPseudoBulk(sce,
  assay = "counts",
  cluster_id = "cell",
  sample_id = "id",
  verbose = FALSE
)
```

## Continuous variable

Consider the continuous variable `Age`. Typical analysis only considers
linear effects using a single regression coefficient, but we also want
to consider the non-linear effects of age. We can peform a [basis
expansion using
splines](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-019-0666-3)
instead use 3 coefficients to model the age effect.

``` r

# Simulate age between 18 and 65
pb$Age <- runif(ncol(pb), 18, 65)

# formula included non-linear effects of Age
# by using a natural spline of degree 3
# This corresponds to using 3 coefficients instead of 1
form <- ~ splines::ns(Age, 3)

# Normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, form, min.count = 5)

# Differential expression analysis within each assay
res.dl <- dreamlet(res.proc, form)

# The spline has degree 3, so there are 3 coefficients
# estimated for Age effects
coefNames(res.dl)
```

    ## [1] "(Intercept)"          "splines::ns(Age, 3)1" "splines::ns(Age, 3)2"
    ## [4] "splines::ns(Age, 3)3"

``` r

# Jointly test effects of the 3 spline components
# The test of the 3 coefficients is performed with an F-statistic
topTable(res.dl, coef = coefNames(res.dl)[2:4], number = 3)
```

    ## DataFrame with 3 rows and 9 columns
    ##         assay                     ID splines..ns.Age..3.1 splines..ns.Age..3.2
    ##   <character>            <character>            <numeric>            <numeric>
    ## 1 CD4 T cells                  GTF3A             0.751933              2.24124
    ## 2 CD4 T cells                   RGS2            -0.591826             -3.72216
    ## 3 CD4 T cells HLA-DRB1_ENSG0000019..            -0.924649             -1.94493
    ##   splines..ns.Age..3.3   AveExpr         F     P.Value adj.P.Val
    ##              <numeric> <numeric> <numeric>   <numeric> <numeric>
    ## 1            -0.418741   8.63391   16.4736 2.82474e-05  0.178952
    ## 2             0.385533   6.97805   14.9336 5.15218e-05  0.178952
    ## 3             1.251838   5.24781   14.7922 5.45632e-05  0.178952

## Ordered categorical

We can also test non-linear effects in the case of categorical variables
with a natural ordering to the categories. Consider time course data
with 4 time points. Each time point is a category and has a natural
ordering from first to last.

We have multiple options to model the time course.

- **Continuous:** Modeling time point as a continuous variable uses a
  single regression coefficient to model the linear effects of the time
  course. This is simple, models the order of the time points, but
  ignores non-linear effects

  Model using `as.numeric(TimePoint)`

- **Categorical:** Including time point as a typical categorical
  variable uses estimated the mean response value for each category. So
  it estimates 4 coefficients. While this can be useful for comparing
  two categories, it ignores the order of the time points.

  Model using `factor(TimePoint)`

- **Ordered categorical:** Here, the trend across ordered time points is
  modled using orthogonal polynomials. The trend is decomposed into
  independent linear, quadratic, etc., effects that can be tested either
  jointly or by themselves.

  Model using:

  ``` r

  ord <- c("time_1", "time_2", "time_3", "time_4")
  ordered(factor(TimePoint), ord)
  ```

Here we simulated 4 time points, and perform differential expression
analysis.

``` r

# Consider data generated across 4 time points
# While there are no time points in the real data
# we can add some for demonstration purposes
pb$TimePoint <- ordered(paste0("time_", rep(1:4, 4)))

# examine the ordering
pb$TimePoint
```

    ##  [1] time_1 time_2 time_3 time_4 time_1 time_2 time_3 time_4 time_1 time_2
    ## [11] time_3 time_4 time_1 time_2 time_3 time_4
    ## Levels: time_1 < time_2 < time_3 < time_4

``` r

# Use formula including time point
form <- ~TimePoint

# Normalize and apply voom/voomWithDreamWeights
res.proc <- processAssays(pb, form, min.count = 5)

# Differential expression analysis within each assay
res.dl <- dreamlet(res.proc, form)

# Examine the coefficient estimated
# for TimePoint it estimates
# linear (i.e. L)
# quadratic (i.e. Q)
# and cubic (i.e. C) effects
coefNames(res.dl)
```

    ## [1] "(Intercept)" "TimePoint.L" "TimePoint.Q" "TimePoint.C"

``` r

# Test only linear effect
topTable(res.dl, coef = "TimePoint.L", number = 3)
```

    ## DataFrame with 3 rows and 9 columns
    ##         assay          ID     logFC   AveExpr         t     P.Value adj.P.Val
    ##   <character> <character> <numeric> <numeric> <numeric>   <numeric> <numeric>
    ## 1 CD4 T cells        DCXR -0.671645   6.52867  -5.42880 4.78666e-05  0.393058
    ## 2 CD4 T cells        GGA2 -0.890112   4.98987  -5.26372 6.69370e-05  0.393058
    ## 3 CD8 T cells        FTH1 -0.759875  14.55380  -4.99199 8.13672e-05  0.393058
    ##           B     z.std
    ##   <numeric> <numeric>
    ## 1  1.686147  -5.42880
    ## 2  0.936906  -5.26372
    ## 3  1.627753  -4.99199

``` r

# Test linear, quadratic and cubic effcts
coefs <- c("TimePoint.L", "TimePoint.Q", "TimePoint.C")
topTable(res.dl, coef = coefs, number = 3)
```

    ## DataFrame with 3 rows and 9 columns
    ##         assay                   ID TimePoint.L TimePoint.Q TimePoint.C
    ##   <character>          <character>   <numeric>   <numeric>   <numeric>
    ## 1 CD8 T cells                 CD52    0.903032   -1.528456   -0.984065
    ## 2 CD8 T cells CCL5_ENSG00000161570    0.626201   -0.986502   -0.645468
    ## 3 CD8 T cells                  CD2    0.887578   -1.152415   -0.195616
    ##     AveExpr         F     P.Value adj.P.Val
    ##   <numeric> <numeric>   <numeric> <numeric>
    ## 1   8.69458   16.0946 1.90868e-05 0.0853298
    ## 2  11.88811   16.0433 1.94982e-05 0.0853298
    ## 3   9.67181   15.7793 2.17777e-05 0.0853298

### Sample filtering

Due to variation in cell and read count for each sample,
[`processAssays()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/processAssays.md)
filters out some sample. This filtering is summarized here:

``` r

details(res.dl)
```

    ##               assay n_retain    formula formDropsTerms n_genes n_errors
    ## 1           B cells       16 ~TimePoint          FALSE    1961        0
    ## 2   CD14+ Monocytes       16 ~TimePoint          FALSE    3087        0
    ## 3       CD4 T cells       16 ~TimePoint          FALSE    5262        0
    ## 4       CD8 T cells       16 ~TimePoint          FALSE    1030        0
    ## 5   Dendritic cells       13 ~TimePoint          FALSE     164        0
    ## 6 FCGR3A+ Monocytes       16 ~TimePoint          FALSE    1160        0
    ## 7    Megakaryocytes       13 ~TimePoint          FALSE     172        0
    ## 8          NK cells       16 ~TimePoint          FALSE    1656        0
    ##   error_initial
    ## 1         FALSE
    ## 2         FALSE
    ## 3         FALSE
    ## 4         FALSE
    ## 5         FALSE
    ## 6         FALSE
    ## 7         FALSE
    ## 8         FALSE

Whle all 16 samples are detained in B cells, only 9 are retained for
megakaryocytes. This can result in a time point being dropped, and so
the polynomial expansion for some cell types can have a lower degree.
The combined results will then have `NA` values for these coefficients.
For example, for `TIMP1` in `Megakaryocytes` above there is not enought
data to fit the cubic term, so `TimePoint.C` is `NA`.

## Session Info

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.30/lib/libopenblasp-r0.3.30.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] muscData_1.24.0             scater_1.38.0              
    ##  [3] scuttle_1.20.0              ExperimentHub_3.0.0        
    ##  [5] AnnotationHub_4.0.0         BiocFileCache_3.0.0        
    ##  [7] dbplyr_2.5.1                muscat_1.24.0              
    ##  [9] dreamlet_1.9.1              SingleCellExperiment_1.32.0
    ## [11] SummarizedExperiment_1.40.0 Biobase_2.70.0             
    ## [13] GenomicRanges_1.62.1        GenomeInfoDb_1.46.2        
    ## [15] Seqinfo_1.0.0               IRanges_2.44.0             
    ## [17] S4Vectors_0.48.0            BiocGenerics_0.56.0        
    ## [19] generics_0.1.4              MatrixGenerics_1.22.0      
    ## [21] matrixStats_1.5.0           variancePartition_1.40.2   
    ## [23] BiocParallel_1.44.0         limma_3.66.0               
    ## [25] ggplot2_4.0.1               BiocStyle_2.38.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fs_1.6.6                  bitops_1.0-9             
    ##   [3] httr_1.4.7                RColorBrewer_1.1-3       
    ##   [5] doParallel_1.0.17         Rgraphviz_2.54.0         
    ##   [7] numDeriv_2016.8-1.1       sctransform_0.4.3        
    ##   [9] tools_4.5.1               backports_1.5.0          
    ##  [11] R6_2.6.1                  metafor_4.8-0            
    ##  [13] mgcv_1.9-4                GetoptLong_1.1.0         
    ##  [15] withr_3.0.2               gridExtra_2.3            
    ##  [17] prettyunits_1.2.0         fdrtool_1.2.18           
    ##  [19] cli_3.6.5                 textshaping_1.0.4        
    ##  [21] sandwich_3.1-1            slam_0.1-55              
    ##  [23] sass_0.4.10               KEGGgraph_1.70.0         
    ##  [25] SQUAREM_2021.1            mvtnorm_1.3-3            
    ##  [27] S7_0.2.1                  blme_1.0-7               
    ##  [29] pkgdown_2.2.0             mixsqp_0.3-54            
    ##  [31] systemfonts_1.3.1         zenith_1.12.0            
    ##  [33] dichromat_2.0-0.1         parallelly_1.46.1        
    ##  [35] invgamma_1.2              RSQLite_2.4.5            
    ##  [37] shape_1.4.6.1             gtools_3.9.5             
    ##  [39] dplyr_1.1.4               Matrix_1.7-4             
    ##  [41] metadat_1.4-0             ggbeeswarm_0.7.3         
    ##  [43] abind_1.4-8               lifecycle_1.0.5          
    ##  [45] yaml_2.3.12               edgeR_4.8.2              
    ##  [47] mathjaxr_2.0-0            gplots_3.3.0             
    ##  [49] SparseArray_1.10.8        grid_4.5.1               
    ##  [51] blob_1.3.0                crayon_1.5.3             
    ##  [53] lattice_0.22-7            beachmat_2.26.0          
    ##  [55] msigdbr_25.1.1            annotate_1.88.0          
    ##  [57] KEGGREST_1.50.0           pillar_1.11.1            
    ##  [59] knitr_1.51                ComplexHeatmap_2.26.0    
    ##  [61] rjson_0.2.23              boot_1.3-32              
    ##  [63] corpcor_1.6.10            future.apply_1.20.1      
    ##  [65] codetools_0.2-20          glue_1.8.0               
    ##  [67] data.table_1.18.0         vctrs_0.7.1              
    ##  [69] png_0.1-8                 Rdpack_2.6.5             
    ##  [71] gtable_0.3.6              assertthat_0.2.1         
    ##  [73] cachem_1.1.0              zigg_0.0.2               
    ##  [75] xfun_0.56                 rbibutils_2.4.1          
    ##  [77] S4Arrays_1.10.1           Rfast_2.1.5.2            
    ##  [79] reformulas_0.4.3.1        iterators_1.0.14         
    ##  [81] statmod_1.5.1             nlme_3.1-168             
    ##  [83] pbkrtest_0.5.5            bit64_4.6.0-1            
    ##  [85] filelock_1.0.3            progress_1.2.3           
    ##  [87] EnvStats_3.1.0            bslib_0.9.0              
    ##  [89] TMB_1.9.19                irlba_2.3.5.1            
    ##  [91] vipor_0.4.7               KernSmooth_2.23-26       
    ##  [93] otel_0.2.0                colorspace_2.1-2         
    ##  [95] rmeta_3.0                 DBI_1.2.3                
    ##  [97] DESeq2_1.50.2             tidyselect_1.2.1         
    ##  [99] bit_4.6.0                 compiler_4.5.1           
    ## [101] curl_7.0.0                httr2_1.2.2              
    ## [103] graph_1.88.1              BiocNeighbors_2.4.0      
    ## [105] desc_1.4.3                DelayedArray_0.36.0      
    ## [107] bookdown_0.46             scales_1.4.0             
    ## [109] caTools_1.18.3            remaCor_0.0.20           
    ## [111] rappdirs_0.3.4            stringr_1.6.0            
    ## [113] digest_0.6.39             minqa_1.2.8              
    ## [115] rmarkdown_2.30            aod_1.3.3                
    ## [117] XVector_0.50.0            RhpcBLASctl_0.23-42      
    ## [119] htmltools_0.5.9           pkgconfig_2.0.3          
    ## [121] lme4_2.0-0                sparseMatrixStats_1.22.0 
    ## [123] lpsymphony_1.38.0         mashr_0.2.79             
    ## [125] fastmap_1.2.0             rlang_1.1.7              
    ## [127] GlobalOptions_0.1.3       htmlwidgets_1.6.4        
    ## [129] UCSC.utils_1.6.1          DelayedMatrixStats_1.32.0
    ## [131] farver_2.1.2              jquerylib_0.1.4          
    ## [133] IHW_1.38.0                zoo_1.8-15               
    ## [135] jsonlite_2.0.0            BiocSingular_1.26.1      
    ## [137] RCurl_1.98-1.17           magrittr_2.0.4           
    ## [139] Rcpp_1.1.1                viridis_0.6.5            
    ## [141] babelgene_22.9            EnrichmentBrowser_2.40.0 
    ## [143] stringi_1.8.7             MASS_7.3-65              
    ## [145] plyr_1.8.9                listenv_0.10.0           
    ## [147] parallel_4.5.1            ggrepel_0.9.6            
    ## [149] Biostrings_2.78.0         splines_4.5.1            
    ## [151] hms_1.1.4                 circlize_0.4.17          
    ## [153] locfit_1.5-9.12           ScaledMatrix_1.18.0      
    ## [155] reshape2_1.4.5            BiocVersion_3.22.0       
    ## [157] XML_3.99-0.20             evaluate_1.0.5           
    ## [159] RcppParallel_5.1.11-1     BiocManager_1.30.27      
    ## [161] nloptr_2.2.1              foreach_1.5.2            
    ## [163] tidyr_1.3.2               purrr_1.2.1              
    ## [165] future_1.69.0             clue_0.3-66              
    ## [167] scattermore_1.2           ashr_2.2-63              
    ## [169] rsvd_1.0.5                broom_1.0.11             
    ## [171] xtable_1.8-4              fANCOVA_0.6-1            
    ## [173] viridisLite_0.4.2         ragg_1.5.0               
    ## [175] truncnorm_1.0-9           tibble_3.3.1             
    ## [177] lmerTest_3.2-0            glmmTMB_1.1.14           
    ## [179] memoise_2.0.1             beeswarm_0.4.0           
    ## [181] AnnotationDbi_1.72.0      cluster_2.1.8.1          
    ## [183] globals_0.18.0            GSEABase_1.72.0
