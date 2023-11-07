
# # TODO
# # Estimateimate sigmaSq by ML
# # fastest for DealyedArray

#' @export
getWeightFromCounts = function(countMatrix){

	count.gene <- rowSums2(countMatrix, useNames=FALSE)
	count.lib <- colSums2(countMatrix, useNames=FALSE)
	ncell <- ncol(countMatrix)
	sclSq <- sum(count.lib^2)

	# add pseudocount
	count.gene <- count.gene + 0.25
	count.lib <- count.lib + 1

	# normalize counts by library size
	# add pseudocount to counts here
	normCounts <- scale(countMatrix + 0.25, 
					scale = count.lib, 
					cente = FALSE)
	# compute variance for each row
	sigmaSq.hat.gene <- rowVars(normCounts, useNames=FALSE)
	sigmaSq.hat.gene[is.na(sigmaSq.hat.gene)] <- 0

	# compute variance
	# vectorize
	v.hat <- 1 / count.gene + (sigmaSq.hat.gene * sclSq) / (ncell^2 *count.gene^2)
	1 / v.hat
}

#' @export
getWeightsForCellType = function(sce, cluster_id, sample_id, CT, weightCap){

	W <- lapply( unique(sce[[sample_id]]), function(ID){

		idx <- sce[[cluster_id]] == CT & sce[[sample_id]] == ID
		countMatrix = counts(sce)[,idx,drop=FALSE]

		getWeightFromCounts( countMatrix )
		})
	W <- do.call(cbind, W)

	# scale each gene by min value
	# so lowest weight is now 1
	W <- W / rowMins( W, useNames=FALSE)

	# set a weight cap at weightCap
	W[W > weightCap] <- weightCap

	# scale each gene to have a mean of 1
	W <- W / rowMeans2(W, useNames=FALSE)

	colnames(W) <- unique(sce[[sample_id]])
	rownames(W) <- rownames(sce)
	W
}

#' @export
getWeightsList = function(sce, cluster_id, sample_id, weightCap = 10){

	if( ! cluster_id %in% colnames(colData(sce)) ){
		msg <- paste0("sample_id entry not found in colData(sce): ", cluster_id)
		stop( msg )
	}
	if( ! sample_id %in% colnames(colData(sce)) ){
		msg <- paste0("sample_id entry not found in colData(sce): ", sample_id)
		stop( msg )
	}

	W.list <- lapply( unique(sce[[cluster_id]]), function(CT){
		getWeightsForCellType( sce, cluster_id, sample_id, CT, weightCap)
	})
	names(W.list)<- unique(sce[[cluster_id]])
 	W.list
}



# #' @export
# trimWeightOutliersGene = function(x, zmax){

# 	# compute z-score
# 	zscore = scale(x)

# 	# extract parameters of transform
# 	# z-score = (x - mu) / s
# 	mu = attr(zscore,"scaled:center")
# 	s = attr(zscore,"scaled:scale")
	
# 	# if x exceeds original value giving z-score of zmax, 
# 	# replace with that orginal value
# 	x[x > zmax * s + mu] = zmax * s + mu

# 	# normalize values to have a mean of 1
# 	x / mean(x)
# }

# #' @export
# trimWeightOutliers = function(W, zmax = 5){

# 	t(apply(W, 1, trimWeightOutliersGene, zmax = zmax))
# }

# # bootstraps cells

countMatrix <- assay(pb, CT)
    	lib.size <- colSums2(countMatrix)
    	names(lib.size) = colnames(countMatrix)
    	# prior.count <- lambda * lib.size/mean(lib.size)
    	df = data.frame(prior.count, n = tab[names(prior.count),CT])
    	prior.count = lambda * df$n
    	prior.count = prior.count / mean(prior.count)
    	names(prior.count) = names(lib.size)


getBootLCPM = function(sce, cluster_id, sample_id, df_pc, ndraws = NULL){
    # interate thu donors, cell types and bootstrap reps
    df_grid = expand.grid(cellType = unique(sce[[cluster_id]]),
                        ID =  unique(sce[[sample_id]]))

    # bootstrap indeces
    idx = sapply( seq(nrow(df_grid)), function(i){

        # filter
        idx = which(df_grid$cellType[i] == sce[[cluster_id]] & df_grid$ID[i] == sce[[sample_id]])

        # bootstrap cells
        if( is.null(ndraws) ){
            idx2 = idx[sample.int(length(idx), length(idx), replace=TRUE)]
        }else{          
            idx2 = idx[sample.int(length(idx), min(length(idx), ndraws), replace=TRUE)]
        }
        idx2
        })
    idx = sort(unlist(idx))

    # pseudobulk of boostrap
    pb <- aggregateToPseudoBulk(sce[,idx],
      assay = "counts",
      cluster_id = cluster_id,
      sample_id = sample_id,
      verbose = FALSE)

    geneExpr = lapply( assayNames(pb), function(CT){

    	df_sub = df_pc %>% filter(cellType == CT)
    	pc = df_sub$prior.count
    	names(pc) = df_sub$ID
    		
    	countMatrix = assay(pb, CT)
    	pcMat = lapply(colnames(countMatrix), function(id)
    				rpois(nrow(countMatrix), pc[id]))
    	pcMat = do.call(cbind, pcMat)
    	colnames(pcMat) = colnames(countMatrix)

        dge = DGEList(counts = countMatrix + pcMat, 
        					lib.size = colSums2(countMatrix))
        # dge = calcNormFactors(dge)
        edgeR::cpm(dge, log=TRUE, prior.count=.25)
        })
    names(geneExpr) = assayNames(pb)

    geneExpr
}

summarizeBootstraps = function(geneExprBoot){
    # interate thu donors, cell types and bootstrap reps
    CT.names = names(geneExprBoot[[1]])
    id.names = colnames(geneExprBoot[[1]][[1]])

    df_var = lapply( CT.names, function(CT){

        df_var = lapply(id.names, function(id){

            # create matrix of boostrap samples for cell type and id
            Y = lapply( seq(length(geneExprBoot)), function(j){
                geneExprBoot[[j]][[CT]][,id,drop=FALSE]
            })
            Y = do.call(cbind, Y)

            # y.mean = rowMeans2(Y, useNames=FALSE)

            # sampling variance of mean from boostraps
            y.var = rowVars(Y, useNames=TRUE) / ncol(Y)

            y.var = data.frame(var = y.var)
            colnames(y.var) = id

            y.var
        })
        as.matrix(do.call(cbind, df_var))
    })
    names(df_var) = CT.names
    df_var
}

library(edgeR)

lambda = 2

CT = "Megakaryocytes"
sceSub = sce[,sce$cell == CT]
colData(sceSub) = droplevels(colData(sceSub))

# options(dplyr.summarise.inform = FALSE)

tab = table(sceSub$id, sceSub$cell)

# prior count *per cell* scaled by variation in library size
lib.size <- colSums2(counts(sceSub))
prior.count <- lambda * lib.size/mean(lib.size)
df_pc = data.frame(ID = sceSub[['id']], 
	cellType = sceSub[['cell']], 
	prior.count = prior.count) %>%
	group_by(cellType, ID) %>%
	summarize(prior.count = sum(prior.count), n=length(ID))


nboots = 103
res = lapply(seq(nboots), function(i) getBootLCPM(sceSub, "cell", "id", df_pc ))


res2 = summarizeBootstraps(res)

V = res2[[CT]]

# n cells
tab = t(with(colData(sce), table(cell, id)))

df = (tab[colnames(V),CT,drop=FALSE])
df = matrix(df, nrow=nrow(V), ncol=length(df), byrow=TRUE)
df[] = nboots

res = limma::squeezeVar( as.numeric(V), as.numeric(df))

par(mfrow=c(1,3))
hist(V)
hist(res$var.post)

plot(as.numeric(V), res$var.post, log="xy")
abline(0, 1, col='red')



A = 
nboots = 103
res = lapply(seq(nboots), function(i) getBootLCPM(sceSub, "cell", "id", df_pc ))


res2 = summarizeBootstraps(res)

V = res2[[CT]]

# n cells
tab = t(with(colData(sce), table(cell, id)))

df = (tab[colnames(V),CT,drop=FALSE])
df = matrix(df, nrow=nrow(V), ncol=length(df), byrow=TRUE)
df[] = nboots

res = limma::squeezeVar( as.numeric(V), as.numeric(df))

par(mfrow=c(1,3))
hist(V)
hist(res$var.post)

plot(as.numeric(V), res$var.post, log="xy")
abline(0, 1, col='red')

A = V / rowMeans(V)



A = V / rowMeans(V)






# n cells
tab = with(colData(sce), table(cell, id))

CT = "B cells"

rv = rowVars(res2[[CT]], useNames=FALSE)
i = which(rv > 0)

df = data.frame(var = res2[[CT]][4,])
df = merge(df, data.frame(n = tab[CT,]), by="row.names")

with(df, plot(n, var))
# x = seq(1, 500)
# lines(x, 1/(lambda*x))

# ratio
max(df$var) / min(df$var)
w = 1/df$var
w = w / mean(w)
hist(w)


df = cellCounts(sce)[,"B cells",drop=FALSE]
V = res2[[1]]
res = limma::squeezeVar( as.numeric(V), 75)

par(mfrow=c(1,3))
hist(V)
hist(res$var.post)

plot(as.numeric(V), res$var.post)
abline(0, 1, col='red')







values = sapply(i, function(j){
	df = data.frame(var = res2[[CT]][j,])
	df = merge(df, data.frame(n = tab[CT,]), by="row.names")
	with(df, cor(n, var, method="sp"))
})


hist(values)






res = lapply(seq(200), function(i) getBootLCPM(sce, "cluster_id", "sample_id"))



res2 = summarizeBootstraps(res)


df = data.frame(n = 2:500)
nboot = 2033
lib.size = 50005
df$var = sapply( df$n, function(n){

	s = sapply(seq(nboot), function(i){
		value = sum(rpois(n, .25)) / (lib.size*n)
		value = value / 1e6
		log(value)
	})
	var(s)
})


with(df, plot(n, var))
# lines(df$n, lambda*df$n / (lib.size * df$n)^2 / 1e12, col="red")
lines(df$n, 1 / (lambda*df$n) , col="red")


# weights proportional to N if 
# Poisson with equal rate
with(df, plot(n, 1/var))
# lines(df$n, lambda*df$n / (lib.size * df$n)^2 / 1e12, col="red")
lines(df$n, lambda*df$n , col="red")



n = 100
lambda = 0.25

n*lambda / lib.size*n

lambda / lib.size



x = rpois(1, lambda*n)
lambda*n / (lib.size * n)^2





    	# average pseudocount of 2
    	pc = rpois(nrow(countMatrix), lambda)

        dge = DGEList(counts = countMatrix + pc,
        					lib.size = colSums(countMatrix))
        # dge = calcNormFactors(dge)



