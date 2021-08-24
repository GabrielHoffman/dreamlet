# Gabriel Hoffman
# August 24, 2021

# Adapted from muscat and scuttle to allow faster access of 
# SingleCellExperiment backed by H5AD file.  
# Uses DelayedMatrixStats for faster computation of summary statistics

.pb = function (x, by, assay, fun, BPPARAM = SerialParam()) 
{
    y <- summarizeAssayByGroup(x, assay.type = assay, ids = (ids <- colData(x)[by]), 
        statistics = fun, BPPARAM = BPPARAM)
    colnames(y) <- y[[by[length(by)]]]
    if (length(by) == 1) 
        return(assay(y))
    if (is.factor(ids <- y[[by[1]]])) 
        ids <- droplevels(ids)
    is <- split(seq_len(ncol(y)), ids)
    ys <- purrr::map(is, ~assay(y)[, .])
    for (i in seq_along(ys)) {
        fill <- setdiff(unique(y[[by[2]]]), colnames(ys[[i]]))
        if (length(fill != 0)) {
            foo <- matrix(0, nrow(x), length(fill))
            colnames(foo) <- fill
            foo <- cbind(ys[[i]], foo)
            o <- paste(sort(unique(y[[by[2]]])))
            ys[[i]] <- foo[, o]
        }
    }
    return(ys)
}



aggregateToPseudoBulk = function (x, assay = NULL, by = c("cluster_id", "sample_id"), 
    fun = c("sum", "mean", "median", "prop.detected", "num.detected"), 
    scale = FALSE, verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)){
    fun <- match.arg(fun)
    if (is.null(assay)) 
        assay <- assayNames(x)[1]
    muscat:::.check_arg_assay(x, assay)
    muscat:::.check_args_aggData(as.list(environment()))
    stopifnot(is(BPPARAM, "BiocParallelParam"))
    for (i in by) if (!is.factor(x[[i]])) 
        x[[i]] <- factor(x[[i]])
    pb <- .pb(x, by, assay, fun, BPPARAM)
    if (scale & length(by) == 2) {
        cs <- if (assay == "counts" && fun == "sum") 
            pb
        else .pb(x, by, "counts", "sum", BPPARAM)
        ls <- lapply(cs, colSums)
        pb <- lapply(seq_along(pb), function(i) pb[[i]]/1e+06 * 
            ls[[i]])
        names(pb) <- names(ls)
    }
    md <- metadata(x)
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







##########################
##########################

#' @importFrom BiocParallel SerialParam 
#' @importFrom SummarizedExperiment SummarizedExperiment
.summarize_assay_by_group <- function(x, ids, subset.row=NULL, subset.col=NULL,
    statistics=c("mean", "sum", "num.detected", "prop.detected", "median"),
    store.number="ncells", threshold=0, BPPARAM=SerialParam()){
    
    new.ids <- .process_ids(x, ids, subset.col)

    sum.out <- .summarize_assay(x, ids=new.ids, subset.row=subset.row,
        statistics=match.arg(statistics, several.ok=TRUE),
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

#' @importFrom BiocParallel SerialParam 
#' @importFrom beachmat rowBlockApply
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

    # out <- rowBlockApply(x, FUN=.summarize_assay_internal, by.group=by.group, 
    #     statistics=statistics, threshold=threshold, BPPARAM=BPPARAM)

    # collected <- do.call(mapply, c(list(FUN=rbind, SIMPLIFY=FALSE, USE.NAMES=FALSE), out))
    # names(collected) <- names(out[[1]])

    freq <- lengths(by.group)

    # if ("mean" %in% statistics) {
    #     collected$mean <- t(t(collected$sum)/freq)
    # }
    # if ("prop.detected" %in% statistics) {
    #     collected$prop.detected <- t(t(collected$num.detected)/freq)
    # }

    # loop through groups

    # Original version uses rowBlockApply() and is slow
    # Use matrixStats and DelayedMatrixStats 
    resCombine = pblapply( by.group, function(idx){

        # subset data by column
        dataSub = x[,idx,drop=FALSE]

        resLst = list()

        # evaluate statistics       
        if( "mean" %in% statistics ){
            resLst[["mean"]] = rowMeans2(dataSub)
        }

        if( "sum" %in% statistics ){
            resLst[["sum"]] = rowSums2(dataSub)
        }

        if( "num.detected" %in% statistics ){
            resLst[["num.detected"]] = rowSums(dataSub > threshold)
        }

        if( "prop.detected" %in% statistics ){          
            resLst[["prop.detected"]] = (ncol(dataSub) - rowCounts(dataSub, value=0)) / ncol(dataSub)
        }

        if( "median" %in% statistics ){
            resLst[["median"]] = rowMedians(dataSub)
        }

        resLst
    }, BPPARAM=BPPARAM)

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

    list(summary=collected, freq=freq)
}

#' @importFrom Matrix rowSums
#' @importFrom DelayedMatrixStats rowMedians 
#' @importClassesFrom Matrix sparseMatrix
#' @importClassesFrom DelayedArray SparseArraySeed
.summarize_assay_internal <- function(x, by.group, statistics, threshold) {
    if (is(x, "SparseArraySeed")) {
        x <- as(x, "sparseMatrix")
    }

    collated <- list()
    
    if ("sum" %in% statistics || "mean" %in% statistics) {
        out <- lapply(by.group, function(i) rowSums(x[,i,drop=FALSE]))
        collated$sum <- .cbind_empty(out, x)
    }

    if ("median" %in% statistics) {
        out <- lapply(by.group, function(i) rowMedians(x[,i,drop=FALSE]))
        out <- .cbind_empty(out, x)
        rownames(out) <- rownames(x)
        collated$median <- out
    }

    if ("num.detected" %in% statistics || "prop.detected" %in% statistics) {
        out <- lapply(by.group, function(i) rowSums(x[,i,drop=FALSE] > threshold))
        collated$num.detected <- .cbind_empty(out, x)
    }

    collated
}

.cbind_empty <- function(out, x) {
    if (length(out)) {
        do.call(cbind, out)
    } else {
        as.matrix(x[,0,drop=FALSE])
    }
}

##########################
##########################

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

#' @export
#' @rdname summarizeAssayByGroup2
setGeneric("summarizeAssayByGroup2", function(x, ...) standardGeneric("summarizeAssayByGroup2"))

#' @export
#' @rdname summarizeAssayByGroup2
#' @importFrom BiocParallel SerialParam
setMethod("summarizeAssayByGroup2", "ANY", .summarize_assay_by_group)

#' @export
#' @rdname summarizeAssayByGroup2
#' @importFrom SummarizedExperiment assay
setMethod("summarizeAssayByGroup2", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .summarize_assay_by_group(assay(x, assay.type), ...)
})







