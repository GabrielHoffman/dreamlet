

#' Forest plot
#'
#' Forest plot
#'
#' @param x result from \code{dreamlet}
#' @param gene gene to show results for
#' @param coef coefficient to test with \code{topTable}
#' @param assays array of assays to plot
#' @param ylim limits for the y axis
#' @param ... other arguments
#'
#' @return Plot showing effect sizes
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
#' # voom-style normalization
#' res.proc = processAssays( pb, ~ group_id)
#' 
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data 
#' res.dl = dreamlet( res.proc, ~ group_id)
#'
#' # show coefficients estimated for each cell type
#' coefNames(res.dl)
#' 
#' # Show estimated log fold change with in each cell type
#' plotForest(res.dl, gene = "ISG20", coef = "group_idstim")
#' 
#' @rdname plotForest-methods
#' @export
setGeneric('plotForest', function(x, gene, coef,...){
	standardGeneric("plotForest")
	})





#' @rdname plotForest-methods
#' @aliases plotForest,dreamletResult-method
#' @export
setMethod("plotForest", signature(x="dreamletResult"),
	function(x, gene, coef, assays=names(x), ylim=NULL){

	# Pass R CMD check
	Assay = logFC = FDR = se = NULL

 	df = topTable(x, coef=coef, number=Inf)
 	df$se = df$logFC / df$t
	df$FDR = p.adjust(df$P.Value)
	df = as.data.frame(df)
	df = df[df$assay %in% assays,]
	df$assay = factor(df$assay, assays)

	ggplot(df[df$ID == gene, ], aes(assay, logFC,  color=-log10(pmax(1e-4,FDR)) )) + 
		geom_point() + geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0.1) + 
		theme_classic() + 
		theme(plot.title = element_text(hjust = 0.5)) + 
		coord_flip(ylim=ylim) + 
		ggtitle(gene) + ylab(bquote(log[2]~fold~change)) + 
		geom_hline(yintercept=0, linetype="dashed") + 
		scale_color_gradient(name = bquote(-log[10]~FDR), low="grey", high="red", limits=c(0, 4)) + 
		xlab('')
})





#' @rdname plotForest-methods
#' @aliases plotForest,dreamlet_mash_result-method
#' @importFrom ashr get_pm get_lfsr get_psd
#' @export
setMethod("plotForest", signature(x="dreamlet_mash_result"),
	function(x, gene, coef, assays=colnames(x$logFC.original), ylim=NULL){

	# Pass R CMD check
	ID = logFC = lFSR = se = NULL

	# extract mashr statistics
	df_logFC = reshape2::melt(get_pm(x$model))
	colnames(df_logFC) = c("Gene", "ID", "logFC")
	df_logFC$key = paste(df_logFC$Gene, df_logFC$ID)

	df_lfsr = reshape2::melt(get_lfsr(x$model))
	colnames(df_lfsr) = c("Gene", "ID", "lFSR")
	df_lfsr$key = paste(df_lfsr$Gene, df_lfsr$ID)

	df_se = reshape2::melt(get_psd(x$model))
	colnames(df_se) = c("Gene", "ID", "se")
	df_se$key = paste(df_se$Gene, df_se$ID)

	# merge mashr statistics
	df = merge(df_logFC, df_lfsr, by='key')
	df = merge(df, df_se, by='key')

	df = df[df$ID %in% assays,]
	df$ID = factor(df$ID, assays)

	# drop NA values
	df = df[!is.na(df$logFC),]

	# make plot
	ggplot(df[df$Gene.x == gene, ], aes(ID, logFC,  color=-log10(pmax(1e-4,lFSR)) )) + 
		geom_point() + 
		geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0.1) + 
		theme_classic() + 
		theme(plot.title = element_text(hjust = 0.5)) + 
		coord_flip(ylim=ylim) + 
		ggtitle(gene) + 
		ylab(bquote(log[2]~fold~change)) + 
		geom_hline(yintercept=0, linetype="dashed") + 
		scale_color_gradient(name = bquote(-log[10]~lFSR), low="grey", high="red", limits=c(0, 4)) + 
		xlab('')
})




