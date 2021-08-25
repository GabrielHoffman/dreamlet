# dreamlet 0.0.5
	August 25, 2021
	- handling of pmetadata by processAssays(), fitVarPart(), and dreamlet()
	- change defaults for bpparam to SerialParam()
	- Created new files for code
	- Create object dreamletResult returned by dreamlet() and used by topTable()
	- more capable dreamletProcessedData object


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