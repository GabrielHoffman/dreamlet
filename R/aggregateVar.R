# Donghoon Lee
# Feb 24, 2023
# modified from aggregateNonCountSignal.R

#' Per-sample variance of single-cell counts
#' 
#' Aggregation function for single-cell log-normalized counts to calculate per-sample variance for dreamlet.  
#' 
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying the assay slot to use as input data. Defaults to the 1st available (\code{assayNames(x)[1]}).
#' @param cluster_id character string specifying which variable to use as cluster id
#' @param sample_id character string specifying which variable to use as sample id
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.var minimum variance for a gene to be considered expressed in a sample
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts for a gene to be
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam}}
#'   object specifying how aggregation should be parallelized.
#' @param verbose logical. Should information on progress be reported?
#'
#' @return a \code{dreamletProcessedData} object 
#' 
#' @details 
#' The \code{dreamlet} workflow can also be applied to model gene expression variance. In this case, a per-sample per-gene variance is calculated across all cells from a given sample and cell type. Here \code{aggregateVar()} performs the roles of \code{aggregateToPseudoBulk()} followed by \code{processAssays()} but using log-normalized count data.
#'
#' For each cell cluster, samples with at least min.cells are retained. Only clusters with at least min.samples retained samples are kept. Features are retained if they have at least min.var in at least min.prop fraction of the samples.
#'  
#' The precision of a measurement is the inverse of its sampling variance. The precision weights are computed as \code{1/sem^2}, where \code{sem = sd / sqrt(n)} and \code{n} is the number of cells.
#' 
#' @examples 
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # Compute variance for each sample and cell cluster
#' pbVar <- aggregateVar(example_sce, 
#'    assay = "counts",    
#'    cluster_id = 'cluster_id', 
#'    sample_id = 'sample_id',
#'    verbose=FALSE)
#' @export
#' @importClassesFrom limma EList
#' @importFrom MatrixGenerics rowVars
aggregateVar = function(sce, assay = NULL, cluster_id = NULL, sample_id = NULL,
                        min.cells = 10, min.var = 0.01,  min.samples = 4, min.prop = 0.4,
                        verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)){

    # get the standard error of the mean for each pseudobulk unit
    pb.sem <- aggregateToPseudoBulk(sce, 
        assay = assay,    
        cluster_id = cluster_id, 
        sample_id = sample_id,
        verbose=FALSE, 
        fun = "sem",
        BPPARAM = BPPARAM,
        checkValues = FALSE)

    # get number of units used for pseudobulk
    pb.number <- aggregateToPseudoBulk(sce, 
        assay = assay,    
        cluster_id = cluster_id, 
        sample_id = sample_id,
        verbose=FALSE, 
        fun = "number",
        BPPARAM = BPPARAM,
        checkValues = FALSE)
    
    # extract metadata shared across assays
    data_constant = droplevels(as.data.frame(colData(pb.sem)))
    
    # variance for each cell type as lists
    resList = lapply(assayNames(pb.sem), function(CT){

        # number of cells used for pseudobulk aggregation
        number = assay(pb.number, CT)

        # variance for each cell type, sem = sd / sqrt(n), var = sd^2
        var = (assay(pb.sem, CT) * sqrt(number))^2

        # precision weights are in inverse variance. so 1/sem^2
        w = 1 / assay(pb.sem, CT)^2

        # create EList with SD and precision weights
        obj = new('EList', list(E = as.matrix(var), weights = as.matrix(w)))

        # inclusion criteria for samples
        include = (number[1,] >= min.cells)
        obj = obj[,include,drop=FALSE]

        # inclusion criteria for genes
        keep = apply(obj$E, 1, function(x){
            sum(x >= min.var) >= min.prop*length(x)
            } )
        obj = obj[keep,,drop=FALSE]

        # if weight is Inf, set to largest finite value per row
        if( any(!is.finite(obj$weights)) ){
            for(i in seq(nrow(obj$weights))){
                if( any(!is.finite(obj$weights[i,])) ){
                    idx = which(!is.finite(obj$weights[i,]))
                    if( length(idx) > 0){
                        obj$weights[i,idx] = max(obj$weights[i, -idx])
                    }
                }
            }
        }

        # if there are too few remaining samples
        if( ncol(obj) < min.samples ){
            return( NULL )
        }

        obj
    })
    names(resList) = assayNames(pb.sem)

    # remove empty assays
    resList = resList[!vapply(resList, is.null, FUN.VALUE=logical(1))]
    
    # return var as dreamletProcessedData object
    new("dreamletProcessedData", resList,
        data = data_constant,
        metadata = metadata(pb.sem)$aggr_means,
        by = c(cluster_id, sample_id))
}
