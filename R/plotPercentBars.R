# Gabriel Hoffman
# Nov 4, 2021


#' Bar plot of variance fractions
#'
#' Bar plot of variance fractions for a subset of genes
#'
#' @param varPart object returned by extractVarPart() or fitExtractVarPartModel()
#' @param col color of bars for each variable
#' @param width specify width of bars
#' 
#' @return Bar plot showing variance fractions for each gene
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
#' # variance partitioning analysis
#' vp = fitVarPart( res.proc, ~ group_id)
#' 
#' # Show variance fractions at the gene-level for each cell type
#' genes = vp$gene[2:4]
#' plotPercentBars(vp[vp$gene %in% genes,])
#' 
#' @export
#' @rdname plotPercentBars-methods
#' @aliases plotPercentBars,vpDF,vpDF-method
#' @importFrom reshape2 melt
#' @importFrom S4Vectors as.data.frame
#' @import ggplot2
setMethod("plotPercentBars", "vpDF",
	function( varPart, col=c(ggColorHue(ncol(varPart)-3), "grey85"), width=NULL){
		
	# convert matrix to tall data.frame
	df = melt( as.data.frame(varPart), id.vars=c("assay", "gene"))

	# pass R CMD check
	gene = value = variable = NA

	fig = ggplot(df, aes(x = gene, y = 100*value, fill = variable)) + 
		geom_bar(stat = "identity", width=width) + theme_bw() + 
		theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank()) + coord_flip() + 
		xlab("") + theme(plot.title=element_text(hjust=0.5)) + facet_wrap(~ assay)

	fig = fig + theme(axis.line = element_line(colour = "transparent"),
		axis.line.x = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(), 
		axis.ticks.y = element_blank(), 
		# legend.position = "bottom",
		# plot.margin = unit(c(0,.3,0,.8), "cm"),
		legend.key = element_blank(),
		panel.spacing.x = unit(.8, "lines")) +
		guides(fill=guide_legend(title=NULL)) +
		scale_fill_manual( values = col) + 
		scale_y_reverse(breaks=seq(0, 100, by=20), label=seq(100, 0, by=-20), expand=c(.00,0)) + 
		ylab("Variance explained (%)")

	fig = fig + theme(text 		= element_text(colour="black"), 
					axis.text 	= element_text(colour="black"),
					legend.text = element_text(colour="black")) 

	fig	
})


