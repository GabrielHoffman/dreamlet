# Gabriel Hoffman
# August 24, 2021

# Adapted from muscat and scuttle to allow faster access of 
# SingleCellExperiment backed by H5AD file.  
# Uses DelayedMatrixStats for faster computation of summary statistics
# This also dramatically reduces memory usage for large datasets
#
# The only change is to .summarize_assay()
# The rest of the code is imported here beccause it is private in muscat, 
# and I didn't want to overwrite summarizeAssayByGroup() in scuttle


#' @importFrom BiocParallel SerialParam
#' @importFrom purrr map
# @importFrom scuttle summarizeAssayByGroup
#' @importFrom SummarizedExperiment assay colData
# @import Matrix
.pb = function (x, by, assay, fun, BPPARAM = SerialParam()) 
{
    y <- summarizeAssayByGroup2(x, assay.type = assay, ids = (ids <- colData(x)[by]), 
        statistics = fun, BPPARAM = BPPARAM)
    colnames(y) <- y[[by[length(by)]]]
    if (length(by) == 1) 
        return(assay(y))
    if (is.factor(ids <- y[[by[1]]])) 
        ids <- droplevels(ids)
    is <- split(seq_len(ncol(y)), ids)
    ys <- purrr::map(is, ~assay(y)[, .])
    for (i in seq_along(ys)) {
        #message(i)
        fill <- setdiff(unique(y[[by[2]]]), colnames(ys[[i]]))
        if (length(fill != 0)) {
            #foo <- matrix(0, nrow(x), length(fill))
            foo <- sparseMatrix(i=1, j=1, x=0, dims=c(nrow(x), length(fill)))
            colnames(foo) <- fill
            foo <- cbind(ys[[i]], foo)
            o <- paste(sort(unique(y[[by[2]]])))
            ys[[i]] <- foo[, o]
        }
    }
    return(ys)
}


# assignInNamespace(".pb", .pb, ns="dreamlet")


#' Aggregation of single-cell to pseudobulk data
#' 
#' Aggregation of single-cell to pseudobulk data.  Adapted from \code{muscat::aggregateData} and has same syntax and results.  But can be much faster in for \code{SingleCellExperiment} backed by H5AD files
#' 
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay character string specifying the assay slot to use as 
#'   input data. Defaults to the 1st available (\code{assayNames(x)[1]}).
# @param by character vector specifying which 
#   \code{colData(x)} columns to summarize by (at most 2!).
#' @param sample_id character string specifying which variable to use as sample id
#' @param cluster_id character string specifying which variable to use as cluster id
#' @param fun a character string.
#'   Specifies the function to use as summary statistic.
#'   Passed to \code{summarizeAssayByGroup2}.
#' @param scale logical. Should pseudo-bulks be scaled
#'   with the effective library size & multiplied by 1M?
#' @param BPPARAM a \code{\link[BiocParallel]{BiocParallelParam}}
#'   object specifying how aggregation should be parallelized.
#' @param verbose logical. Should information on progress be reported?
#' 
#' @return a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
# \itemize{
# \item{If \code{length(by) == 2}, each sheet (\code{assay}) contains 
#   pseudobulks for each of \code{by[1]}, e.g., for each cluster when 
#   \code{by = "cluster_id"}. Rows correspond to genes, columns to 
#   \code{by[2]}, e.g., samples when \code{by = "sample_id"}}.
# \item{If \code{length(by) == 1}, the returned SCE will contain only 
#   a single \code{assay} with rows = genes and colums = \code{by}.}}
#' 
#'  Aggregation parameters (\code{assay, by, fun, scaled}) are stored in  \code{metadata()$agg_pars}, where \code{by = c(cluster_id, sample_id)}.  The number of cells that were aggregated are accessible in \code{int_colData()$n_cells}.
#' 
#' @examples 
#' 
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce, 
#'    assay = "counts",    
#'    cluster_id = 'cluster_id', 
#'    sample_id = 'sample_id',
#'    verbose=FALSE)
#'
#' # pseudobulk data from each cell type
#' # is stored as its own assay 
#' pb
#' 
#' # aggregate by cluster only, 
#' # collapsing all samples into the same pseudobulk
#' pb2 <- aggregateToPseudoBulk(example_sce, cluster_id = "cluster_id", verbose=FALSE)
#' 
#' pb2  
#' 
#' @author Gabriel Hoffman, Helena L Crowell & Mark D Robinson
#'
#' @details 
#' Adapted from \code{muscat::aggregateData} and has simular syntax and same results.  This is much faster for \code{SingleCellExperiment} backed by H5AD files using \code{DelayedMatrix} because this summarizes counts using \code{\link[DelayedMatrixStats]{DelayedMatrixStats}}.  But this function also includes optmizations for \code{sparseMatrix} used by \code{\link[Seurat]{Seurat}} by using \code{sparseMatrixStats}.
#' 
#' @references 
#' Crowell, HL, Soneson, C, Germain, P-L, Calini, D, 
#' Collin, L, Raposo, C, Malhotra, D & Robinson, MD: 
#' On the discovery of population-specific state transitions from 
#' multi-sample multi-condition single-cell RNA sequencing data. 
#' \emph{bioRxiv} \strong{713412} (2018). 
#' doi: \url{https://doi.org/10.1101/713412}
#' 
#' @importFrom Matrix colSums
#' @importFrom purrr map
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims<- int_colData<-
#' @importFrom SummarizedExperiment rowData colData colData<-
#' @export
aggregateToPseudoBulk = function (x, assay = NULL, sample_id = NULL, cluster_id = NULL,
    fun = c("sum", "mean", "median", "prop.detected", "num.detected"), 
    scale = FALSE, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)){
    
    fun <- match.arg(fun)

    # specify the 'by' parameter using cluster_id and sample_id.
    by = c(cluster_id, sample_id)

    if( is.null(by) ){
        stop("Must specify at least one (and usually both): sample_id, cluster_id ")
    }

    if (is.null(assay)) 
        assay <- assayNames(x)[1]
    .check_arg_assay(x, assay)
    .check_args_aggData(as.list(environment()))
    stopifnot(is(BPPARAM, "BiocParallelParam"))
    for (i in by) if (!is.factor(x[[i]])) 
        x[[i]] <- factor(x[[i]])

    # store metdata
    md <- metadata(x)

    # set to empty so not passed to downstream functions
    metadata(x) = list() 
    reducedDims(x) = list()
        
    suppressPackageStartupMessages({ 
        pb <- .pb(x, by, assay, fun, BPPARAM)
    })
    if (scale & length(by) == 2) {
        cs <- if (assay == "counts" && fun == "sum") 
            pb
        else suppressPackageStartupMessages({ 
            .pb(x, by, "counts", "sum", BPPARAM)})
        ls <- lapply(cs, colSums)
        pb <- lapply(seq_along(pb), function(i) pb[[i]]/1e+06 * 
            ls[[i]])
        names(pb) <- names(ls)
    }
    md$agg_pars <- list(assay = assay, by = by, fun = fun, scale = scale)
    pb <- SingleCellExperiment(pb, rowData = rowData(x), metadata = md)
    cd <- data.frame(colData(x)[, by])
    for (i in names(cd)) if (is.factor(cd[[i]])) 
        cd[[i]] <- droplevels(cd[[i]])
    ns <- table(cd)
    if (length(by) == 2) {
        ns <- asplit(ns, 2)
        ns <- purrr::map(ns, ~c(unclass(.)))
    }
    else ns <- c(unclass(ns))
    int_colData(pb)$n_cells <- ns
    if (length(by) == 2) {
        cd <- colData(x)
        ids <- colnames(pb)
        counts <- vapply(ids, function(u) {
            m <- as.logical(match(cd[, by[2]], u, nomatch = 0))
            vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
        }, numeric(ncol(colData(x))))
        cd_keep <- apply(counts, 1, function(u) all(u == 1))
        cd_keep <- setdiff(names(which(cd_keep)), by) 
        if (length(cd_keep) != 0) {
            m <- match(ids, cd[, by[2]], nomatch = 0)
            cd <- cd[m, cd_keep, drop = FALSE]
            rownames(cd) <- ids
            colData(pb) <- cd
        }
    }
    return(pb)
}




#' @importFrom SummarizedExperiment assayNames
.check_arg_assay = function (x, y){
    stopifnot(is.character(y), length(y) == 1)

    if( ! y %in% assayNames(x)){
        stop("Assay name not found: ", y)
    }
    if (sum(assayNames(x) == y) > 1){
        stop("Argument 'assay' was matched to multiple times.\n ", 
            " Please assure that the input SCE has unique 'assayNames'.")
    }

    # check that assay is integers
    #-----------------------------
    # get nonzero values from first 3 samples 
    values = as.matrix(assay(x,y)[,1:3])
    values = values[values!=0]

    if( any(values < 0) ){
        txt = paste0("Assay '", y, "' stores negatives values instead of counts.\n\tThis is not allowed for creating pseudobulk.\n\tThis dataset contains assays: ", paste(assayNames(x), collapse=', '))
        stop(txt)
    }

    # compute fraction that are already integers
    integerFrac = sum(values == floor(values)) / length(values)

    if( integerFrac < 1 ){
        txt = paste0("Assay '", y, "' stores continuous values instead of counts.\n\tMake sure this is the intended assay.\n\tThis dataset contains assays: ", paste(assayNames(x), collapse=', '))
        warning(txt, immediate.=TRUE)
    }
}


#' @importFrom SummarizedExperiment colData
.check_args_aggData = function (u) 
{
    stopifnot(is.character(u$by), length(u$by) <= 2, u$by %in% 
        colnames(colData(u$x)))
    stopifnot(is.logical(u$scale), length(u$scale) == 1)
    if (u$scale & (!u$assay %in% c("cpm", "CPM") | u$fun != "sum")) 
        stop("Option 'scale = TRUE' only valid for", " 'assay = \"cpm/CPM\"' and 'fun = \"sum\"'.")
}




##########################
##########################

#' @importFrom BiocParallel SerialParam 
#' @importFrom SummarizedExperiment SummarizedExperiment
.summarize_assay_by_group <- function(x, ids, subset.row=NULL, subset.col=NULL,
    statistics=c("mean", "sum", "num.detected", "prop.detected", "median"),
    store.number="ncells", threshold=0, BPPARAM=SerialParam()){
    
    statistics = match.arg(statistics, several.ok=TRUE)

    new.ids <- .process_ids(x, ids, subset.col)

    sum.out <- .summarize_assay(x, ids=new.ids, subset.row=subset.row,
        statistics=statistics,
        threshold=threshold, BPPARAM=BPPARAM)

    mat.out <- sum.out$summary
    mapping <- match(colnames(mat.out[[1]]), as.character(new.ids))
    coldata <- .create_coldata(original.ids=ids, mapping=mapping, 
        freq=sum.out$freq, store.number=store.number)

    # Sync'ing the names as coldata is the source of truth.
    for (i in seq_along(mat.out)) {
        colnames(mat.out[[i]]) <- rownames(coldata)
    }

    SummarizedExperiment(mat.out, colData=coldata)
}

#' @importFrom BiocParallel SerialParam bplapply
#' @useDynLib dreamlet
#' @import Rcpp RcppEigen RcppProgress
#' @importFrom DelayedArray colsum
#' @importFrom methods as
.summarize_assay <- function(x, ids, statistics, threshold=0, subset.row=NULL, BPPARAM=SerialParam()) {

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    lost <- is.na(ids)
    ids <- ids[!lost]
    if (any(lost)) {
        x <- x[,!lost,drop=FALSE]
    }

    # Drop unused levels, as the subsequent mapping step to preserve the type of 'ids'
    # in .create_coldata doesn't make sense (as there is no mapping to a concrete observation).
    by.group <- split(seq_along(ids), ids, drop=TRUE)

    # frequency of each cluster type
    freq <- lengths(by.group)

    # method for aggregation depends on datatype    
    #  I used RcppEigen to speed up rowSums for sparseMatrix and matrix
    if( (length(statistics) == 1) & (statistics[1] == "sum") & is(x, "sparseMatrix") ){
        # countsMatrix = rowSums_by_chunk_sparse(x, by.group, verbose=BPPARAM$progressbar)
        # rownames(countsMatrix) = rownames(x)
        # colnames(countsMatrix) = names(by.group)

        countsMatrix = colsum_beachmat(x, ids)    

        collected = list(sum = as(countsMatrix, "sparseMatrix"))
    }else if( (length(statistics) == 1) & (statistics[1] == "sum") & is(x, "matrix") ){
        # countsMatrix = rowSums_by_chunk(x, by.group, verbose=BPPARAM$progressbar)
        # rownames(countsMatrix) = rownames(x)
        # colnames(countsMatrix) = names(by.group)

        countsMatrix = colsum_beachmat(x, ids)          
            
        collected = list(sum = as(countsMatrix, "sparseMatrix"))
    }else if( (length(statistics) == 1) & (statistics[1] == "sum") & is(x, "DelayedMatrix") ){

         # countsMatrix <- colsum_fast(x, ids, BPPARAM=BPPARAM)
        verbose = BPPARAM$progressbar
        BPPARAM$progressbar = FALSE
        countsMatrix <- colsum2(x, ids, BPPARAM=BPPARAM, verbose=BPPARAM$progressbar)        

        collected = list(sum = as(countsMatrix, "sparseMatrix"))
    }else{    
        collected = .pb_summary(x, by.group, statistics, threshold, BPPARAM)
    }

    list(summary=collected, freq=freq)
}


# method for aggregation depends on datatype    
#' @importFrom DelayedMatrixStats rowMeans2 rowSums2 rowCounts rowMedians
#' @importFrom sparseMatrixStats rowMeans2 rowSums2 rowCounts rowMedians
#' @importFrom MatrixGenerics rowMeans2 rowSums2 rowCounts rowMedians
.pb_summary = function(x, by.group, statistics, threshold, BPPARAM){

    if( is(x, "DelayedMatrix") ){
        
        # Original version uses rowBlockApply() and is slow
        # Use DelayedMatrixStats dramatically increases speed
        resCombine = bplapply( by.group, function(idx, data){

            # when run in parallel, each thread loads packages.  So suppress
            suppressPackageStartupMessages({ 

            resLst = list()

            # evaluate statistics       
            if( "mean" %in% statistics ){
                resLst[["mean"]] = DelayedMatrixStats::rowMeans2(data, cols=idx)
            }

            if( "sum" %in% statistics ){
                resLst[["sum"]] = DelayedMatrixStats::rowSums2(data, cols=idx)
            }

            if( "num.detected" %in% statistics ){
                resLst[["num.detected"]] = DelayedMatrixStats::rowSums2(data[,idx,drop=FALSE] > threshold)
            }

            if( "prop.detected" %in% statistics ){          
                resLst[["prop.detected"]] = (length(idx) - DelayedMatrixStats::rowCounts(data, value=0, cols=idx)) / length(idx)
            }

            if( "median" %in% statistics ){
                resLst[["median"]] = DelayedMatrixStats::rowMedians(data, cols=idx)
            }
            })

            resLst
        }, data=x, BPPARAM=BPPARAM)

    }else if( is(x, "sparseMatrix") ){
        # Avoid repeated garbage collection steps in bclapply()
        resCombine = lapply( by.group, function(idx){

            resLst = list()

            # evaluate statistics       
            if( "mean" %in% statistics ){
                resLst[["mean"]] = sparseMatrixStats::rowMeans2(x, cols=idx)
            }

            if( "sum" %in% statistics ){
                resLst[["sum"]] = sparseMatrixStats::rowSums2(x, cols=idx)
            }

            if( "num.detected" %in% statistics ){
                resLst[["num.detected"]] = sparseMatrixStats::rowSums2(x[,idx,drop=FALSE] > threshold)
            }

            if( "prop.detected" %in% statistics ){          
                resLst[["prop.detected"]] = (length(idx) - sparseMatrixStats::rowCounts(x, value=0, cols=idx)) / length(idx)
            }

            if( "median" %in% statistics ){
                resLst[["median"]] = sparseMatrixStats::rowMedians(x, cols=idx)
            }

            resLst
        })
          
    }else{
        resCombine = lapply( by.group, function(idx){

            resLst = list()

            # evaluate statistics       
            if( "mean" %in% statistics ){
                resLst[["mean"]] = rowMeans(x[,idx,drop=FALSE])
            }

            if( "sum" %in% statistics ){
                resLst[["sum"]] = rowSums2(x[,idx,drop=FALSE])
            }

            if( "num.detected" %in% statistics ){
                resLst[["num.detected"]] = rowSums2(x[,idx,drop=FALSE] > threshold)
            }

            if( "prop.detected" %in% statistics ){          
                resLst[["prop.detected"]] = (length(idx) - rowCounts(x[,idx,drop=FALSE], value=0)) / length(idx)
            }

            if( "median" %in% statistics ){
                resLst[["median"]] = rowMedians(x[,idx,drop=FALSE])
            }

            resLst
        })
    }

    # Create list of merged values for each statistic
    collected = lapply(statistics, function(stat){
        tmp = lapply(resCombine, function(res){
            res[[stat]]
        })
        tmpMat = do.call(cbind, tmp)
        rownames(tmpMat) = rownames(x)
        tmpMat
        })
    names(collected) = statistics
    
    collected
}




# #' @importFrom Matrix rowSums
# #' @importFrom DelayedMatrixStats rowMedians 
# #' @importClassesFrom Matrix sparseMatrix
# #' @importClassesFrom DelayedArray SparseArraySeed
# .summarize_assay_internal <- function(x, by.group, statistics, threshold) {
#     if (is(x, "SparseArraySeed")) {
#         x <- as(x, "sparseMatrix")
#     }

#     collated <- list()
    
#     if ("sum" %in% statistics || "mean" %in% statistics) {
#         out <- lapply(by.group, function(i) rowSums(x[,i,drop=FALSE]))
#         collated$sum <- .cbind_empty(out, x)
#     }

#     if ("median" %in% statistics) {
#         out <- lapply(by.group, function(i) rowMedians(x[,i,drop=FALSE]))
#         out <- .cbind_empty(out, x)
#         rownames(out) <- rownames(x)
#         collated$median <- out
#     }

#     if ("num.detected" %in% statistics || "prop.detected" %in% statistics) {
#         out <- lapply(by.group, function(i) rowSums(x[,i,drop=FALSE] > threshold))
#         collated$num.detected <- .cbind_empty(out, x)
#     }

#     collated
# }

# .cbind_empty <- function(out, x) {
#     if (length(out)) {
#         do.call(cbind, out)
#     } else {
#         as.matrix(x[,0,drop=FALSE])
#     }
# }

##########################
##########################

# scuttle
.subset2index = function (subset, target, byrow = TRUE) 
{
    if (is.na(byrow)) {
        dummy <- seq_along(target)
        names(dummy) <- names(target)
    }
    else if (byrow) {
        dummy <- seq_len(nrow(target))
        names(dummy) <- rownames(target)
    }
    else {
        dummy <- seq_len(ncol(target))
        names(dummy) <- colnames(target)
    }
    if (!is.null(subset)) {
        subset <- dummy[subset]
        if (any(is.na(subset))) {
            stop("invalid subset indices specified")
        }
    }
    else {
        subset <- dummy
    }
    unname(subset)
}


#' @importFrom S4Vectors selfmatch
.df_to_factor <- function(ids) {
    o <- order(ids)
    x <- selfmatch(ids[o,,drop=FALSE]) 
    x[o] <- x
    x[Reduce("|", lapply(ids, is.na))] <- NA_integer_
    x
}

#' @importFrom methods is
.has_multi_ids <- function(ids) is(ids, "DataFrame")

.process_ids <- function(x, ids, subset.col) {    
    if (.has_multi_ids(ids)) {
        ids <- .df_to_factor(ids)
    } 
    if (ncol(x)!=length(ids)) {
        stop("length of 'ids' and 'ncol(x)' are not equal")
    }
    if (!is.null(subset.col)) {
        ids[!seq_along(ids) %in% .subset2index(subset.col, x, byrow=FALSE)] <- NA_integer_
    }
    ids
}

#' @importFrom S4Vectors DataFrame
.create_coldata <- function(original.ids, mapping, freq, store.number) {
    if (.has_multi_ids(original.ids)) {
        coldata <- original.ids[mapping,,drop=FALSE]
        rownames(coldata) <- NULL
    } else {
        coldata <- DataFrame(ids=original.ids[mapping])
        rownames(coldata) <- coldata$ids
    }

    if (!is.null(store.number)) {
        coldata[[store.number]] <- unname(freq)
    }
    coldata
}

##########################
##########################

# Adapted form scuttle::summarizeAssayByGroup

# @export
# @rdname summarizeAssayByGroup2
setGeneric("summarizeAssayByGroup2", function(x, ...) standardGeneric("summarizeAssayByGroup2"))

# @export
# @rdname summarizeAssayByGroup2
# @importFrom BiocParallel SerialParam
setMethod("summarizeAssayByGroup2", "ANY", .summarize_assay_by_group)

# @export
# @rdname summarizeAssayByGroup2
#' @importFrom SummarizedExperiment assay
setMethod("summarizeAssayByGroup2", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .summarize_assay_by_group(assay(x, assay.type), ...)
})





