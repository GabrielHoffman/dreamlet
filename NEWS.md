# dreamlet 1.1.16
  - Feb 8, 2024
  - fix bug in `pbWeights()`

# dreamlet 1.1.15
  - Feb 5, 2024
  - Fix bug in call to `eBayes()`
   - in `processAssays()` pass argument `scaledByLib` to `voomWithDreamWeights()`

# dreamlet 1.1.14
  - Jan 29, 2024
  - fix bug in `pbWeights()`
   - smaller pseudo variance
   - limit to only expressed genes by adding `getExprGeneNames()`

# dreamlet 1.1.13
  - Jan 25, 2024
  - improve error reporting in `seeErrors()` and documentation
  - update `outlier()` to compute z-scores.  How returns `data.frame()`
  - add `outlierByAssay()` and `plotPCA()`

# dreamlet 1.1.12
  - Jan 16, 2024
  - `compositePosteriorTest()` allows `exclude` set to be `NULL`

# dreamlet 1.1.11
  - Jan 16, 2024
  - add `meta_analysis()`

# dreamlet 1.1.10
  - Jan 10, 2024
  - `stackAssays()` now includeds `metadata()$aggr_means` correctly
  - add `compositePosteriorTest()`

# dreamlet 1.1.9
  - Jan 3, 2024
  - use `get_metadata_aggr_means()` to extract `aggr_means` when SCE is produced by cbind'ing

# dreamlet 1.1.8
  - Dec 18, 2023
  - fix issue when no genes pass cutoffs
  - fix issue with `aggr_means` in `aggregateToPseudoBulk()`
  - fix bug when `rdf` is low for all genes

# dreamlet 1.1.7
  - Dec 10, 2023
  - add `plotBeeswarm()` 
  - add `rowWeightedVarsMatrix()`
  - bug fixes
  - add `isFullRank()` check in `dreamlet()`
  - handle exceptions in `run_mash()`

# dreamlet 1.1.6
  - Dec 5, 2023
  - `pbWeights()` add argument `maxRatio`

# dreamlet 1.1.5
  - Nov 28, 2023
  - fix edge cases in `pbWeights()` 
  - cell weights is not default in `dreamlet()`

# dreamlet 1.1.4
  - Nov 27, 2023
  - add `pbWeights()` to compute precision weights for pseudobulk counts
  - extend `extractData()` and include it in vignette

# dreamlet 1.1.3
  - Nov 22, 2023
  - add `stackAssays()`
  - add `diffVar()`
  - fix `getVarFromCounts)` so `zeta` is a mean, not a sum

# dreamlet 1.1.2
  - Nov 13, 2023
  - `computeLogCPM()` uses `augmentPriorCount()`

# dreamlet 1.0.1 / 1.1.1
  - bug fix in Bioc 3.18 and devel

# dreamlet 0.99.33
  - `computeLogCPM()` now returns `matrix` instead of `sparseMatrix`

# dreamlet 0.99.32
  - Oct 18, 2023
  - update `variancePartition` version dependency
  - add `getWeightsList()`

# dreamlet 0.99.31
  - Oct 12, 2023
  - fix error reporting in `processAssays()`

# dreamlet 0.99.30
  - Sept 22, 2023
  - Fix issue with precision weights from cell number
  - Include `setAutoBlockSize()` update within `aggregateToPseudoBulk()`

# dreamlet 0.99.28
  - Sept 5, 2023
  - Update error handling for `processAssays()` and `fitVarPart()`

# dreamlet 0.99.26
  - August 18, 2023
  - Update error handling and documentation

# dreamlet 0.99.23
  - August 8, 2023
  - run `styler::style_pkg()`

# dreamlet 0.99.22
  - August 8, 2023
  - `dreamletCompareClusters()` now allows cell-level covariates in response to https://github.com/GabrielHoffman/dreamlet/issues/11
  - Fix code for Bioconductor submission

# dreamlet 0.99.21
  - July 17, 2023
  - Improve functionality and documentation of `dreamlet::residuals()`

# dreamlet 0.99.20
  - June 29, 2023
  - in `processAssays()` use  `voomWithDreamWeights(..., span="auto")` to estimate the lowess tuning parameter 
  - rare error in `merge_metadata()` when a cell type is not observed for all donors.

# dreamlet 0.99.19
  - June 28, 2023
  - in `dreamlet()` fix issue when contrasts are specified and formula includes variable from `metadata()`

# dreamlet 0.99.18
  - June 27, 2023
  - add `assays` argument to `buildClusterTreeFromPB()`

# dreamlet 0.99.14
  - June 16, 2023
  - bug fix in `processAssays()` when assays is dropped

# dreamlet 0.99.13
  - May 31, 2023
  - improved error reporting
  - Compatibility with variancePartition v2.0.5 (renamed 1.31.1)

# dreamlet 0.99.12
  - May 24, 2023
  - required `zellkonverter (>= 1.10.1)` to avoid issues with previous version
    - issue solved by https://github.com/theislab/zellkonverter/blob/b56718d113327020c024e188d9ac67ea57eaf35d/R/AnnData2SCE.R#L351

# dreamlet 0.99.11
  - May 12, 2023
  - Compatibility with variancePartition v2.0.1

# dreamlet 0.99.10
  - April 24, 2023
  - Fix issue in raised in Bioconductor submission: https://github.com/Bioconductor/Contributions/issues/2955#issuecomment-1498070980
  - Compatibility with variancePartition v2.0.0

# dreamlet 0.99.6
  - March 29, 2023
  - fix `topTable()` for `dreamletResult` in the case where one or more cells didn't estimate the coefficient of interest

# dreamlet 0.99.3
  - March 23, 2023
  - reduce compiler time
  - add `computeNormCounts()` and `computeLogCPM()` 

# dreamlet 0.99.1
  - March 20, 2023
  - fix Biocondcutor submission based on https://github.com/Bioconductor/Contributions/issues/2955#issuecomment-1476037237

# dreamlet 0.99.0
  - March 15, 2023
  - Biocondcutor submission

# dreamlet 0.0.64
 - March 7, 2023
 - fix `topTable()` to deal with multiple `coef` as array
 - add vignette about non-linear effects
 - Add `colData<-` for `dreamletProcessedData`

# dreamlet 0.0.63
 - March 3, 2023
 - small bug fixes in `topTable()` and `plotForest()`

# dreamlet 0.0.62
 - Feb 27, 2023
 - `aggregateToPseudoBulk()` stores mean of cell-level covariates in `metadata(pb)$aggr_means`
  - Use these covariates in `processAssays()`, `dreamlet()`, `fitVarPart()`
  - extend to `aggregateNonCountSignal()`
 
# dreamlet 0.0.61
 - Jan 25, 2023
 - add `plotProjection()`
 - add `outlier()`
 - update `plotForest()`

# dreamlet 0.0.60
 - Jan 18, 2023
 - allow custom ordering of assays in plots

# dreamlet 0.0.59
 - Jan 12, 2023
 - update `plotVolcano()` to allow `scales="free_y"`
   - warning when p-values are zero.

# dreamlet 0.0.58
 - Jan 9, 2023
 - Rewrite`aggregateNonCountSignal()` to include filters
 - Depend only on CRAN and BioC packages

# dreamlet 0.0.57
 - Jan 4, 2023
 - Include precision weights in `aggregateNonCountSignal()`

# dreamlet 0.0.56
 - Jan 3, 2023
 - Enable processing of non-count data with `aggregateNonCountSignal()`

# dreamlet 0.0.55
 - Dec 7, 2022
 - in `plotGeneHeatmap()` drop empty genes

# dreamlet 0.0.54
 - Nov 30, 2022
 - add `buildClusterTreeFromPB()`

# dreamlet 0.0.53
 - Nov 18, 2022
 - fix bug in `topTable()`

# dreamlet 0.0.52
  - Nov 18, 2022
  - add `as.dreamletResult()`
  - update `variancePartition` dependency and source

# dreamlet 0.0.51
  - Nov 10, 2022
  - in `processAssays()` and `processOneAssay()`, add argument `min.prop` indicating the minimum proportion of retained samples with non-zero counts

# dreamlet 0.0.50
  - Nov 2, 2022
  - fix bug in `zenith_gsa()` for few gene sets

# dreamlet 0.0.49
  - Oct 17, 2022
  - add `computeCellCounts()`

# dreamlet 0.0.48
  - Sept 14, 2022
  - add `transpose` argument to `plotGeneHeatmap()`
  - and `alpha` arugment to `plotVoom()`
  - update y axis and outlines in `plotVarPart()`

# dreamlet 0.0.47
  - update filtering of covariates, especially for when many samples are dropped 

# dreamlet 0.0.46
  - add `totalCPM` column to output of `cellTypeSpecificity()` to use for filtering.  Functions `dreamlet::plotHeatmap()` `plotViolin()` and `plotPercentBars()` now ignore this column 

# dreamlet 0.0.45
  - update documentation

# dreamlet 0.0.44
  - add `plotGeneHeatmap()`
  - add argument `assays` to `plotVarPart()`
  - add `extractData()`

# dreamlet 0.0.43
  - faster `aggregateToPseudoBulk()` by speeding up check in `.check_arg_assay()`
  - more flexibility for `tabToMatrix()`

# dreamlet 0.0.42
  - fix misc issues with plotting and order of facets
  - fix issue with redundant variables in small sample sizes

# dreamlet 0.0.41
  - fix issue with `topTable()` when all random effects are dropped

# dreamlet 0.0.40
  - Compatibility with zenith package submitted to Bioconductor.

# dreamlet 0.0.39
  - fix issue with `aggregateToPseudoBulk()` when summarizing for just 1 sample

# dreamlet 0.0.38
  - add `getTreat()` for `dreamlet()` result

# dreamlet 0.0.37
  - `droplevels` for `colData` in `processAssays()`

# dreamlet 0.0.36
  - fixes in `processAssays()` to detect issues with SCE

# dreamlet 0.0.35
  - remove old code and dependencies
  - change website theme

# dreamlet 0.0.34
 - bug fix for `aggregateToPseudoBulk()` with sample ordering

# dreamlet 0.0.33
 - bug fix for `dreamletCompareClusters()`

# dreamlet 0.0.32
 - add `colsum2()` using beachmat code.

# dreamlet 0.0.31
 - reduce memory usage in `aggregateToPseudoBulk()` by fixing `aggregateByColnames()

# dreamlet 0.0.30
 - update `run_mash()` to combine results across coefs

# dreamlet 0.0.29
 - fix bug in `dreamlet::colsum_fast()` used in pseudobulk

# dreamlet 0.0.28
 - add `da_to_sparseMatrix()`

# dreamlet 0.0.27
 - `aggregateToPseudoBulk()` for `DelayedArray` now uses `colsum_fast()`
   - this is faster then the previous version for `DelayedArray`

# dreamlet 0.0.26
  - update `dreamletCompareClusters()`:
    - now compatable with `plotZenithResults()`
    - include flag `errorsAsWarnings`.  If `TRUE` warns and returns NULL. 

# dreamlet 0.0.25
  - change return value for `dreamletCompareClusters()` to be compatible with `zenith_gsa()`
  - fix usage of `formula` in `dreamletCompareClusters()`

# dreamlet 0.0.24
  - additional speed up `aggregateToPseudoBulk()` when a Seurat object is used
    - uses RcppEigen sparse matrix iterator

# dreamlet 0.0.23
  - dramatic speed up `aggregateToPseudoBulk()` when a Seurat object is used
    - uses RcppEigen

# dreamlet 0.0.22
  - Speed up `aggregateToPseudoBulk()` when a Seurat object is used

# dreamlet 0.0.21
  - Add `collapse=TRUE` to `dreamletCompareClusters()`

# dreamlet 0.0.20
  - Fix bug in `dreamletCompareClusters()`

# dreamlet 0.0.19
  - Fix bug in `dreamletCompareClusters()`

# dreamlet 0.0.18
 - add `min.samples` to `processAssays()`, `processOneAssay()`
 - add `dreamletCompareClusters()` and `run_mash()`
  - Fix bug in `dreamletCompareClusters()`
  - updated `mashr` dependency


# dreamlet 0.0.17
 - add `run_mash()`
   * add `zenith_gsa()`, `plotVolcano()`, `plotForest()` for results
 - fix bug in `cellTypeSpecificity()` for genes with zero reads across all cell types
 - order of arguments in `plotForest()` and `zenith_gsa()` changed for consistancy
 - expand vignettes
 - bug fix for `removeConstantTerms()` when excluded variable string (i.e. tissue) is also a substring of other variables (i.e. tissueStatus)

# dreamlet 0.0.16
- add `residuals()` for `dreamlet()` result
- add `dreamletPairs()`
- fix bug in `removeConstantTerms()` with multiple constant terms
- improve usability of `cellTypeSpecificity()` by adding `plotPercentBars()` and `plotViolin()` compatability
- fix bug in `topTable()` when `coef` is not estimated
- add argument `assays` to `dreamlet()`, `fitVarPart()`, and `processAssays()`

# dreamlet 0.0.15
- `processOneAssay()` weights by number of cells
- require `variancePartition >= 1.25.1` to handle weights in `voomWithDreamWeights()`
- fix bug in `topTable()`
- add `plotPercentBars()` for class `vpDF`
- misc bug fixes
- improve documentation


# dreamlet 0.0.14
- move count ratio code to crumblr package
- use `applyQualityWeights()`

# dreamlet 0.0.13
- Oct 25, 2021
- add `ilr_composition_test.R`

# dreamlet 0.0.12
- Oct 15, 2021
- update print for `dreamletResult` using `coefNames()`
- small bug fix
- fix bugs in `regModel()`

# dreamlet 0.0.11
- Oct 6, 2021
- `removeConstantTerms()` now doesn't drop terms solely because of NA's
	- this means that other functions can gracefully warn the user about NA's

# dreamlet 0.0.10
- Oct 5, 2021
- suppress package startup messages in `aggregateToPseudoBulk()`
- bug fix in `removeConstantTerms()`

- Sept 30, 2021
- call to `zenith_gsa()` adds argument `inter.gene.cor` and `progressbar`
- fix output to `cellTypeCompositionVarPart()` and `cellTypeCompositionTest()`
- fix issue with `topTable()` where FDR was evaluated on only a subset of genes

# dreamlet 0.0.8
- Sept 28, 2021
- update docs, logos, TODO
- update `dreamlet()` to handle linear contrasts
- `removeConstantTerms()` now drops categorical variables with only a max of one example per category
- Cleaner code for cell composition test

# dreamlet 0.0.7
- Sept 2, 2021
- add `cellTypeCompositionTest()`
	- handle random effects


# dreamlet 0.0.6
- Sept 1, 2021
- enforce package version requirements

# dreamlet 0.0.5
- August 25, 2021
- handling of pmetadata by processAssays(), fitVarPart(), and dreamlet()
- change defaults for bpparam to SerialParam()
- Created new files for code
- Create object dreamletResult returned by dreamlet() and used by topTable()
- more capable dreamletProcessedData object
- add details() functions for dreamletResult, dreamletProcessedData and vpDF
- More information about when terms are dropped from formulas
- type definition to zenith_gsa


# dreamlet 0.0.4
- add aggregateToPseudoBulk() for faster access to SingleCellExperiment() backed by H5AD

# dreamlet 0.0.3
- add argument to `processAssays()` to include extra meta-data	
- add subseting with assay() for dreamletProcessedData
- adapt plotVoom(), plotVolcano(), plotVarPart() to be more flexiable
- fitVarPart() returns DataFrame


# dreamlet 0.0.2
- check failed regression

# dreamlet 0.0.1
- Initial version