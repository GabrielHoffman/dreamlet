


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
#' @param assays which assays to plot
#' @param ... arguments passed to \code{facet_wrap()}. Useful for specifying \code{scales = "free_y"}
#'
#' @return Volcano plot for each cell type
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
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data 
#' res.dl = dreamlet( res.proc, ~ group_id)
#' 
#' # show coefficients estimated for each cell type
#' coefNames(res.dl)
#' 
#' # volcano plot for each cell type
#' plotVolcano(res.dl, coef="group_idstim")
#'
#' # volcano plot for first two cell types
#' plotVolcano(res.dl[1:2], coef="group_idstim")
#'
#' @export
#' @docType methods
#' @rdname plotVolcano-methods
setGeneric("plotVolcano", 
  function(x, coef, nGenes=5, size=12, minp=1.0e-310, cutoff=0.05, ncol=3,...){

  standardGeneric("plotVolcano")
})



#' @importFrom data.table as.data.table
#' @importFrom ggrepel geom_text_repel
#' @rdname plotVolcano-methods
#' @aliases plotVolcano,list,list-method
setMethod("plotVolcano", "list",
  function(x, coef, nGenes=5, size=12, minp=1.0e-310, cutoff=0.05, ncol=3, assays = names(x),...){

  # intersect preserving order from assays
  assays = intersect(assays, names(x))
  if( length(assays) == 0) stop("No valid assays selected")
    
  df_combine = topTable(x, coef=coef, number=Inf)
  df_combine = as.data.table(df_combine)
  idx = df_combine$assay %in% assays
  df_combine = df_combine[idx,,drop=FALSE]

  # Pass R CMD check
  .SD = logFC = P.Value = isSignif = Gene = ID = NULL

  xmax = max(abs(df_combine$logFC))
  ymax = -log10(min(df_combine$P.Value))

  # check for 0 p-values
  if( ! is.finite(ymax) ){
    nzero = length(df_combine$P.Value==0)
    txt = paste0("There are ", nzero, " features with p-value of 0. Plotting will be affected")
    warning(txt)
  }

  df_combine$isSignif = c("no","yes")[(df_combine$adj.P.Val < cutoff)+1]
  df_combine$P.Value = pmax(minp, df_combine$P.Value )

  # sort facets by original sorting of assays
  df_combine$assay = factor(df_combine$assay, assays)
  # order by assay, then p-value
  ord = order(df_combine$assay, df_combine$P.Value)
  df_combine = df_combine[ord,,drop=FALSE]

  # top significant genes in each cell type
  df2 = df_combine[,head(.SD, nGenes), by="assay",drop=FALSE]

  # reverse order to plot significant points last
  fig = ggplot(df_combine, aes(logFC, -log10(P.Value), color=isSignif)) + 
    geom_point() + 
    theme_classic(size) + 
    theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    xlab(bquote(log[2]~fold~change)) + 
    ylab(bquote(-log[10]~P)) + 
    scale_color_manual(values=c("grey", "darkred")) + 
    geom_text_repel(data=df2, aes(logFC, -log10(P.Value), label=ID), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5) +
    facet_wrap(~assay, ncol=ncol,...) 

  # If scales is not specified, standardize y-axis 
  # and don't expand
  argsList = match.call(expand.dots = FALSE)$`...`
  ags = names(argsList)
  if( length(ags) == 0 || !("scales" %in% ags) ){
    fig = fig + 
    scale_y_continuous(limits=c(0, ymax*1.02), expand=c(0,0))
  }  

  fig
})



#' @importFrom data.table data.table
#' @importFrom ggrepel geom_text_repel
#' @rdname plotVolcano-methods
#' @aliases plotVolcano,MArrayLM,MArrayLM-method
setMethod("plotVolcano", "MArrayLM",
  function(x, coef, nGenes=5, size=12, minp=1.0e-310, cutoff=0.05, ncol=3,...){

  tab = topTable(x, coef=coef, number=Inf)
  df_combine = data.table(ID = rownames(tab), tab)

  xmax = max(abs(df_combine$logFC))
  ymax = -log10(min(df_combine$P.Value))

  # Pass R CMD check
  .SD = logFC = P.Value = isSignif = ID = NULL

  # check for 0 p-values
  if( ! is.finite(ymax) ){
    nzero = length(df_combine$P.Value==0)
    txt = paste0("There are ", nzero, " features with p-value of 0. Plotting will be affected")
    warning(txt)
  }

  df_combine$isSignif = c("no","yes")[(df_combine$adj.P.Val < cutoff)+1]
  df_combine$P.Value = pmax(minp, df_combine$P.Value )

  # top significant genes in each cell type
  df2 = df_combine[,head(.SD, nGenes)]

  # reverse order to plot significant points last
  ggplot(df_combine[seq(nrow(df_combine), 1),,drop=FALSE], aes(logFC, -log10(P.Value), color=isSignif)) + 
    geom_point() + 
    theme_classic(size) + 
    theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    xlab(bquote(log[2]~fold~change)) + 
    ylab(bquote(-log[10]~P)) + 
    scale_color_manual(values=c("grey", "darkred")) + 
    scale_y_continuous(limits=c(0, ymax*1.02), expand=c(0,0)) + 
    geom_text_repel(data=df2, aes(logFC, -log10(P.Value), label=ID), segment.size=.5,  segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5)
})





#' @importFrom data.table data.table
#' @importFrom ggrepel geom_text_repel
#' @importFrom reshape2 melt
#' @rdname plotVolcano-methods
#' @aliases plotVolcano,dreamlet_mash_result,dreamlet_mash_result-method
setMethod("plotVolcano", "dreamlet_mash_result",
  function(x, coef, nGenes=5, size=12, minp=1.0e-16, cutoff=0.05, ncol=3, assays = colnames(x$logFC.original),...){

  # intersect preserving order from assays
  assays = intersect(assays, colnames(x$logFC.original))
  if( length(assays) == 0) stop("No valid assays selected")

  df_logFC = reshape2::melt(get_pm(x$model))
  colnames(df_logFC) = c("Gene", "ID", "logFC")
  df_logFC$key = paste(df_logFC$Gene, df_logFC$ID)

  df_lfsr = reshape2::melt(get_lfsr(x$model))
  colnames(df_lfsr) = c("Gene", "ID", "lFSR")
  df_lfsr$key = paste(df_lfsr$Gene, df_lfsr$ID)

  df = merge(df_logFC, df_lfsr, by='key')
  df = data.table(df[!is.na(df$logFC),])

  # sort by lFSR
  df = df[order(df$lFSR),,drop=FALSE]

  # sort facets by original sorting of assays
  df$ID.x = factor(df$ID.x, colnames(x$logFC.original))

  # Pass R CMD check
  .SD = logFC = P.Value = isSignif = Gene.x = lFSR = NULL

  df$isSignif = c("no","yes")[(df$lFSR < cutoff)+1]
  df$lFSR = pmax(minp, df$lFSR )

  # filter based on assays
  df = df[df$ID.x %in% assays,,drop=FALSE]
  df$ID.x = factor(df$ID.x, assays)

  xmax = max(abs(df$logFC))
  ymax = -log10(min(df$lFSR))

  # top significant genes in each cell type
  df2 = df[,head(.SD, nGenes), by="ID.x"]

  # order by assay, then lFSR
  ord = order(df$ID.x, df$lFSR)
  df = df[ord,,drop=FALSE]

  fig = ggplot(df, aes(logFC, -log10(lFSR), color=isSignif)) + 
    geom_point() + 
    theme_classic(size) + 
    theme(aspect.ratio=1, legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    xlab(bquote(log[2]~fold~change)) + 
    ylab(bquote(-log[10]~local~False~Sign~Rate~(mashr))) + 
    scale_color_manual(values=c("grey", "darkred")) + 
    geom_text_repel(data=df2, aes(logFC, -log10(lFSR), label=Gene.x), segment.size=.5, segment.color="black", color="black", force=1, nudge_x=.005, nudge_y=.5) +
    facet_wrap(~ID.x, ncol=ncol,...) 

  # If scales is not specified, standardize y-axis 
  # and don't expand
  argsList = match.call(expand.dots = FALSE)$`...`
  ags = names(argsList)
  if( length(ags) == 0 || !("scales" %in% ags) ){
    fig = fig + 
    scale_y_continuous(limits=c(0, ymax*1.02), expand=c(0,0))
  }  

  fig
})








