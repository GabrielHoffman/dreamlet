
#' Forest plot
#'
#' Forest plot
#'
#' @param x result from \code{dreamlet}
#' @param coef coefficient to test with \code{topTable}
#' @param gene gene to show results for
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
#' plotForest(res.dl, coef = "group_idstim", gene = "ISG20")
#' 
#' @import variancePartition limma
#' @import ggplot2
#' @export
plotForest = function(x, coef, gene){

	# Pass R CMD check
	Assay = logFC = FDR = se = NULL


	df = lapply( names(x), function(assay){
		tab = topTable(x[[assay]], coef=coef, number=Inf)
		data.frame(Assay = assay, Gene = rownames(tab), tab, se = with(tab, logFC / t))
		})
	df = do.call(rbind, df)
	df$FDR = p.adjust(df$P.Value)

	ggplot(df[df$Gene == gene, ], aes(Assay, logFC,  color=-log10(pmax(1e-4,FDR)) )) + geom_point() + geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0.1) + theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + coord_flip() + ggtitle(gene) + ylab(bquote(log[2]~fold~change)) + geom_hline(yintercept=0, linetype="dashed") + scale_color_gradient(name = bquote(-log[10]~FDR), low="grey", high="red", limits=c(0, 4))
}





