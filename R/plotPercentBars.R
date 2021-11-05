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
#' @export
#' @rdname plotVarPart-methods
#' @aliases plotVarPart,vpDF,vpDF-method
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
		legend.key = element_blank(),
		plot.margin = unit(c(0,.3,0,.8), "cm")) +
		guides(fill=guide_legend(title=NULL)) +
		scale_fill_manual( values = col) + 
		scale_y_reverse(breaks=seq(0, 100, by=20), label=seq(100, 0, by=-20), expand=c(0,0.03)) + 
		ylab("Variance explained (%)")

	fig = fig + theme(text 		= element_text(colour="black"), 
					axis.text 	= element_text(colour="black"),
					legend.text = element_text(colour="black")) 

	fig	
})


