
setClass("vpDF", contains="DFrame", slots=c(df_details = "data.frame"))

#' Bar plot of variance fractions
#'
#' Bar plot of variance fractions for a subset of genes
#'
#' @param x \code{vpDF} object returned by \code{fitVarPart()}
#' @param col color of bars for each variable
#' @param genes name of genes to plot
#' @param width specify width of bars
#' @param ncol number of columns in the plot
#' @param ... other arguments
#' 
#' @return Bar plot showing variance fractions for each gene
#' 
#' @examples
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
#' plotPercentBars(vp, genes = vp$gene[2:4], ncol=2)
#' 
#' @export
#' @rdname plotPercentBars-methods
#' @aliases plotPercentBars,vpDF,vpDF-method
#' @importFrom reshape2 melt
#' @import ggplot2
setMethod("plotPercentBars", "vpDF",
	function( x, col=c(ggColorHue(ncol(x)-3), "grey85"), genes=unique(x$gene), width=NULL, ncol = 3,...){
	
	# get assays when it is not defined in generic function		
	args <- list(...)
	if( 'assays' %in% names(args) ){
		assays = args$assays
	}else{
		assays = assayNames(x)
	}

	# subset based on assays
	x = x[x$assay %in% unique(assays),]	

	# subset based on specified genes
	x = x[x$gene %in% unique(genes),]	

	# convert matrix to tall data.frame
	df = melt( as.data.frame(x), id.vars=c("assay", "gene"))

	df$gene = factor(df$gene, rev(genes))

	# pass R CMD check
	gene = value = variable = NA

	fig = ggplot(df, aes(x = gene, y = 100*value, fill = variable)) + 
		geom_bar(stat = "identity", width=width) + theme_bw() + 
		theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank()) + coord_flip() + 
		xlab("") + theme(plot.title=element_text(hjust=0.5)) + facet_wrap(~ assay, ncol=ncol)

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





#' @export
#' @rdname plotPercentBars-methods
#' @aliases plotPercentBars,cellSpecificityValues,cellSpecificityValues-method
#' @importFrom reshape2 melt
#' @import ggplot2
setMethod("plotPercentBars", "cellSpecificityValues",
	function( x, col=ggColorHue(ncol(x)), genes=rownames(x), width=NULL,...){
		
	gene = unique(genes)
	idx = match(genes, rownames(x))
	idx = idx[!is.na(idx)]

	if( length(idx) == 0 ){
		stop("No genes left after filtering")
	}

	if( length(idx) != length(genes) ){
		txt = paste(length(genes) - length(idx), "genes were not found")
		warning(txt)
	}

	# omit column totalCPM, if it exists
	i = which(colnames(x) == "totalCPM")
	if( length(i) > 0) x= x[,-1]

	df = data.frame(x[idx,], check.names=FALSE)
	df$gene = rownames(df)
	df_melt = melt(df, id.vars="gene")

	df_melt$gene = factor(df_melt$gene, rev(gene))

	# pass R CMD check
	gene = value = variable = NA

	ggplot(df_melt, aes(gene, value, fill=variable)) + 
	  geom_bar(stat="identity") + 
	  theme_classic() + 
	  theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + 
	  coord_flip() + 
	  scale_fill_discrete(name = "Cell type") + 
	  xlab("Gene") + 
	  ylab("Fraction of gene expression") + 
	  scale_y_continuous(limits=c(0, 1 + 1e-14), expand=c(0,0)) + 
	 ggtitle("Cell type specificity scores")
})







