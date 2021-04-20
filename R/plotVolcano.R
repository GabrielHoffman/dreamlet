
#' Volcano plot
#' 
#' Volcano plot
#' 
#' @param df data.frame produced by \code{topTable}
#' @param showGenes array of gene names corresponding rownames of \code{df} to highlight
#' @param size text size
#' @param minp minimum p-value to show on the y-axis
#' @param cutoff adj.P.Val cutoff to distinguish significant from non-significant genes
#' @param xlim limits on the x-axis
#' @param ymax maximum value on the y-axis
#' 
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @export
plotVolcano = function(df, showGenes = NULL, size=15, minp=1e-310, cutoff=0.05, xlim=NULL, ymax=NULL ){

  # PASS R CMD check
  logFC = P.Value = isSignif = Symbol = NULL

  # indicate significant genes
  df$isSignif = c("no","yes")[(df$adj.P.Val < cutoff)+1]
  df$P.Value = pmax(minp, df$P.Value )

  if( is.null(ymax) ) ymax = -log10(min(df$P.Value))

  if( is.null(xlim) ) xlim = c(-max(abs(df$logFC)), max(abs(df$logFC)))

  fig = ggplot(df, aes( logFC, -log10(P.Value), color=isSignif)) + 
    geom_point() + 
    theme_bw(size) + 
    theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    xlab(bquote(log[2]~fold~change)) + 
    ylab(bquote(-log[10]~P)) + 
    scale_color_manual(values=c("grey", "darkred")) + 
    scale_y_continuous(limits=c(0, ymax*1.02), expand=c(0,0)) + 
    xlim(xlim)

  if( !is.null(showGenes) ){
    # top significant genes
    df2 = df[rownames(df) %in% showGenes,]
    df2$Symbol = rownames(df2)

    fig = fig + geom_text_repel(data=df2, aes(logFC, -log10(P.Value), label=Symbol), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5)
  }

  fig
}

#' Volcano plot for each cell type
#' 
#' Volcano plot for each cell type
#' 
#' @param x result from \code{dreamlet}
#' @param coef coefficient to test with \code{topTable}
#' @param nGenes number of genes to highlight in each volcano plot
#' @param size text size
#' @param minp minimum p-value to show on the y-axis
#' @param cutoff adj.P.Val cutoff to distinguish significant from non-significant genes
#' @param ncol number of columns in the plot
#' 
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
plotVolcano_grid = function(x, coef, nGenes=5, size=12, minp=1e-310, cutoff=0.05, ncol=3){

  dfList = lapply( x, function(fit){

    topTable(fit, coef=coef, number=Inf)
    })
  names(dfList) = names(x)
  df_combine = do.call(rbind, dfList)

  xmax = max(abs(df_combine$logFC))
  ymax = -log10(min(df_combine$P.Value))

  figList = lapply( names(dfList), function(cellType){

    df = dfList[[cellType]]

    plotVolcano(df, showGenes=rownames(df)[seq_len(nGenes)], size=size, minp=minp, cutoff=cutoff, 
      xlim = c(-xmax, xmax), ymax=ymax ) + ggtitle(cellType)
    })

  plot_grid(plotlist = figList, ncol=ncol, align="hv", axis='tblr')
}







