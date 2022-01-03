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
Sept 1, 2021
- enforce package version requirements

# dreamlet 0.0.5
August 25, 2021
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