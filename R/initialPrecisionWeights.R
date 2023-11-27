

#' @importFrom dplyr tibble bind_rows
getVarFromCounts <- function(countMatrix, lib.size, prior.count = .25){

    stopifnot( ncol(countMatrix) == length(lib.size))

    countMatrix <- countMatrix + prior.count

    # pseudobulk
    count.gene <- rowSums2(countMatrix, useNames=FALSE)

    # normalize counts by library size
    # add pseudocount to counts here
    normCounts <- scale(countMatrix, 
                    scale = lib.size + 1, 
                    center = FALSE)

    # compute variance for each row
    sigSq.hat <- rowVars(normCounts, useNames=FALSE)
    sigSq.hat[is.na(sigSq.hat)] <- 0

    # return values to compute variance later
    tibble(Gene = rownames(countMatrix), 
            count.gene = count.gene ,
            sigSq.hat = sigSq.hat, 
            zeta = mean(lib.size^2), 
            ncell = ncol(countMatrix))
}


getVarForCellType <- function(sce, sample_id, cluster_id, CT, prior.count){

    cellType <- ID <- NULL

    if( ! "counts" %in% assayNames(sce) ){
        stop("SCE does not contain assay: counts")
    }

    idx <- which(sce[[cluster_id]] == CT)
    lib.size <- colSums2(counts(sce), cols=idx, useNames=TRUE)
    
    # Scale prior count so that an observed count of 0,
    # gives zero variance across samples
    # Add small value to each cell, so that across n_i cells
    # it the augment sum to a mean of prior.count
    df_pc <- data.frame(ID = sce[[sample_id]][idx], 
        cellType = sce[[cluster_id]][idx], 
        prior.count = prior.count * lib.size/mean(lib.size)) %>%
        group_by(cellType, ID) %>%
        summarize(n=length(ID), prior.count = mean(prior.count) / n, .groups="drop_last")

    # get variance estimates for each ID and gene
    df <- lapply( unique(sce[[sample_id]]), function(ID){

        idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == ID
        countMatrix <- counts(sce)[,idx,drop=FALSE]

        pc <- df_pc$prior.count[df_pc$ID == ID]
        
        res <- getVarFromCounts( countMatrix, 
                lib.size = lib.size[colnames(countMatrix)], 
                prior.count = pc)
        res$ID <- ID
        res
        })
    bind_rows(df)
}

#' @importFrom limma squeezeVar
#' @importFrom Matrix sparseMatrix
#' @importFrom dplyr mutate
getVarList <- function(sce, sample_id, cluster_id, shrink = TRUE, prior.count = 0.5){

    Gene = ID = count.gene = ncell = zeta = sigSq.hat = NULL

    if( ! cluster_id %in% colnames(colData(sce)) ){
        msg <- paste0("sample_id entry not found in colData(sce): ", cluster_id)
        stop( msg )
    }
    if( ! sample_id %in% colnames(colData(sce)) ){
        msg <- paste0("sample_id entry not found in colData(sce): ", sample_id)
        stop( msg )
    }

    # Compute variance for each observation for each cell type
    var.list <- lapply( unique(sce[[cluster_id]]), function(CT){

        df <- getVarForCellType( sce, sample_id, cluster_id, CT, prior.count) %>%
                mutate(Gene = factor(Gene, rownames(sce)),
                        ID = factor(ID))

        if( shrink ){      
            # shrink sample variances     
            res <- squeezeVar( df$sigSq.hat, df$ncell-1, robust=FALSE)
            df$sigSq.hat <- res$var.post
        }

        # delta approximation of variance
        df <- df %>%
            mutate( vhat = 1 / count.gene * (1 + ncell*sigSq.hat*zeta / count.gene)) 

        mat <- sparseMatrix(df$Gene, df$ID, 
            x = df$vhat, 
            dims = c(nlevels(df$Gene), nlevels(df$ID)),
            dimnames = list(levels(df$Gene), levels(df$ID)))
        as.matrix(mat)
    })
    names(var.list) <- unique(sce[[cluster_id]])

    var.list
}

#' Compute observation weights for pseudobulk
#' 
#' Compute observation weights for pseudobulk using the delta method to approximate the variance of the log2 counts per million considering variation in the number of cells and gene expression variance across cells within each sample.
#' 
#' @param sce \code{SingleCellExperiment} of where \code{counts(sce)} stores the raw count data at the single cell level
#' @param sample_id character string specifying which variable to use as sample id
#' @param cluster_id character string specifying which variable to use as cluster id
#' @param shrink Defaults to \code{TRUE}. Use empirical Bayes variance shrinkage from \code{limma} to shrink estimates of expression variance across cells within each sample
#' @param prior.count Defaults to \code{0.5}. Count added to each observation at the pseudobulk level.  This is scaled but the number of cells before added to the cell level
#' @param quantileOffset Defaults to \code{0.1}. When computing the precision from the variance, regularize the reciprocal by adding a small value to the denominator. For a gene with variances stored in the array \code{x}, add \code{quantile(x, quantileOffset)} before taking the reciprocal.
#' @param h5adBlockSizes set the automatic block size block size (in bytes) for DelayedArray to read an H5AD file.  Larger values use more memory but are faster.
#'
#' @examples
#' library(muscat)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce,
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#' 
#' # Create precision weights for pseudobulk
#' weightsList <- pbWeights(example_sce, 
#'     cluster_id = "cluster_id",
#'     sample_id = "sample_id")
#' 
#' @importFrom stats quantile
#' @importFrom DelayedArray getAutoBlockSize setAutoBlockSize
#' @export
pbWeights <- function(sce, sample_id, cluster_id, shrink = TRUE, prior.count = 0.5, quantileOffset = 0.1, h5adBlockSizes = 1e9){

    # update block size for reading h5ad file from disk
    tmp <- getAutoBlockSize()
    suppressMessages(setAutoBlockSize(h5adBlockSizes))
    on.exit(suppressMessages(setAutoBlockSize(tmp)))

    # compute variances
    var.lst <- getVarList(sce, sample_id, cluster_id, shrink, prior.count)

    # regularize reciprocal with quantile offset
    W.list <- lapply(var.lst, function(x){
        1 / ( x + quantile(x, quantileOffset))
    })

    W.list
}




