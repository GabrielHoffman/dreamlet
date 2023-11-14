

#' @importFrom dplyr tibble bind_rows
getVarFromCounts = function(countMatrix, lib.size, prior.count = .25){

    stopifnot( ncol(countMatrix) == length(lib.size))

    countMatrix = countMatrix + prior.count

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
            zeta = sum(lib.size^2), 
            ncell = ncol(countMatrix))
}


getVarForCellType = function(sce, cluster_id, sample_id, CT, prior.count){

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
#' @export
getVarList = function(sce, cluster_id, sample_id, shrink, prior.count = 0.5){

    if( ! cluster_id %in% colnames(colData(sce)) ){
        msg <- paste0("sample_id entry not found in colData(sce): ", cluster_id)
        stop( msg )
    }
    if( ! sample_id %in% colnames(colData(sce)) ){
        msg <- paste0("sample_id entry not found in colData(sce): ", sample_id)
        stop( msg )
    }

    # for each cell type
    var.list <- lapply( unique(sce[[cluster_id]]), function(CT){

        df <- getVarForCellType( sce, cluster_id, sample_id, CT, prior.count) %>%
                mutate(Gene = factor(Gene, rownames(sce)),
                        ID = factor(ID))

        if( shrink ){      
            # shrink sample variances     
            res <- squeezeVar( df$sigSq.hat, df$ncell-1, robust=FALSE)
            # plot(df$sigSq.hat, res$var.post, main=CT, log="xy")
            # abline(0, 1, col="red")
            # browser()
            df$sigSq.hat <- res$var.post
        }

        # delta approximation of variance
        df = df %>%
            mutate( vhat = 1 / count.gene * (1 + sigSq.hat*zeta / (count.gene/ncell)))

        mat <- sparseMatrix(df$Gene, df$ID, 
            x = df$vhat, 
            dims = c(nlevels(df$Gene), nlevels(df$ID)),
            dimnames = list(levels(df$Gene), levels(df$ID)))
        as.matrix(mat)
    })
    names(var.list) <- unique(sce[[cluster_id]])
    var.list
}





