


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
#' @import ggplot2 ggrepel
#' @export
#' @docType methods
#' @rdname plotVolcano-methods
setGeneric("plotVolcano", 
  function(x, coef, nGenes=5, size=12, minp=1.0e-310, cutoff=0.05, ncol=3){

  standardGeneric("plotVolcano")
})



#' @rdname plotVolcano-methods
#' @aliases plotVolcano,list,list-method
setMethod("plotVolcano", "list",
  function(x, coef, nGenes=5, size=12, minp=1.0e-310, cutoff=0.05, ncol=3){

  dfList = lapply( names(x), function(id){

    tab = topTable(x[[id]], coef=coef, number=Inf)
    data.table(ID = id, Gene = rownames(tab), tab)
    })
  names(dfList) = names(x)
  df_combine = do.call(rbind, dfList)

  # Pass R CMD check
  .SD = logFC = P.Value = isSignif = Gene = NULL

  xmax = max(abs(df_combine$logFC))
  ymax = -log10(min(df_combine$P.Value))

  df_combine$isSignif = c("no","yes")[(df_combine$adj.P.Val < cutoff)+1]
  df_combine$P.Value = pmax(minp, df_combine$P.Value )

  # top significant genes in each cell type
  df2 = df_combine[,head(.SD, nGenes), by="ID"]

  # reverse order to plot significant points last
  ggplot(df_combine[seq(nrow(df_combine), 1)], aes(logFC, -log10(P.Value), color=isSignif)) + 
    geom_point() + 
    theme_bw(size) + 
    theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    xlab(bquote(log[2]~fold~change)) + 
    ylab(bquote(-log[10]~P)) + 
    scale_color_manual(values=c("grey", "darkred")) + 
    scale_y_continuous(limits=c(0, ymax*1.02), expand=c(0,0)) + 
    geom_text_repel(data=df2, aes(logFC, -log10(P.Value), label=Gene), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5) +
    facet_wrap(~ID, ncol=ncol) 
})




#' @rdname plotVolcano-methods
#' @aliases plotVolcano,MArrayLM,MArrayLM-method
setMethod("plotVolcano", "MArrayLM",
  function(x, coef, nGenes=5, size=12, minp=1.0e-310, cutoff=0.05, ncol=3){

  tab = topTable(x, coef=coef, number=Inf)
  df_combine = data.table(Gene = rownames(tab), tab)

  xmax = max(abs(df_combine$logFC))
  ymax = -log10(min(df_combine$P.Value))


  # Pass R CMD check
  .SD = logFC = P.Value = isSignif = Gene = NULL

  df_combine$isSignif = c("no","yes")[(df_combine$adj.P.Val < cutoff)+1]
  df_combine$P.Value = pmax(minp, df_combine$P.Value )

  # top significant genes in each cell type
  df2 = df_combine[,head(.SD, nGenes)]

  # reverse order to plot significant points last
  ggplot(df_combine[seq(nrow(df_combine), 1)], aes(logFC, -log10(P.Value), color=isSignif)) + 
    geom_point() + 
    theme_bw(size) + 
    theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    xlab(bquote(log[2]~fold~change)) + 
    ylab(bquote(-log[10]~P)) + 
    scale_color_manual(values=c("grey", "darkred")) + 
    scale_y_continuous(limits=c(0, ymax*1.02), expand=c(0,0)) + 
    geom_text_repel(data=df2, aes(logFC, -log10(P.Value), label=Gene), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5)
})











