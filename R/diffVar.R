#' Class dreamletResult
#'
#' Class \code{dreamletResult} stores results produced by \code{dreamlet()} to give a standard interface for downstream analysis
#'
#' @name dreamletResult-class
#' @rdname dreamletResult-class
#' @exportClass dreamletResult
#' @return none
setClass("dreamletResult", contains = "list", slots = c(df_details = "data.frame", errors = "list", error.initial = "list"))



setGeneric("diffVar", variancePartition::diffVar)

#' Test differential variance
#'
#' Test the association between a covariate of interest and the response's deviation from expectation.
#'
#' @param fit model fit from \code{dream()}
#' @param method transform the residuals using absolute deviation ("AD") or squared deviation ("SQ").
#' @param scale scale each observation by "leverage", or no scaling ("none")
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other parameters passed to \code{dream()}
#'
#' @details
#' This method performs a test of differential variance between two subsets of the data, in a way that generalizes to multiple categories, continuous variables and metrics of spread beyond variance.  For the two category test, this method is simular to Levene's test.  This model was adapted from Phipson, et al (2014), extended to linear mixed models, and adapted to be compatible with \code{variancePartition::dream()} and \code{dreamlet::dreamlet()}.
#'
#' This method is composed of multiple steps where 1) a typical linear (mixed) model is fit with \code{dreamlet()}, 2) residuals are computed and transformed based on an absolute value or squaring transform, 3) a second regression is performed with \code{dreamlet()} to test if a variable is associated with increased deviation from expectation.  Both regression take advantage of the \code{dreamlet()} linear (mixed) modelling framework followed by empirical Bayes shrinkage that extends the \code{limma::voom()} framework.
#'
#' Note that \code{diffVar()} takes the results of the first regression as a parameter to use as a starting point.
#'
#' @references{
#'   \insertRef{phipson2014diffvar}{variancePartition}
#' }
#' @seealso \code{variancePartition::diffVar()}
#' @examples
#' library(muscat)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce,
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#'
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~ group_id)
#'
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data
#' res.dl <- dreamlet(res.proc, ~ group_id)
#'
#' # Differential variance analysis
#' # result is a dreamlet fit
#' res.dvar <- diffVar( res.dl )
#' 
#' # Examine results
#' res.dvar
#'
#' # Examine details for each assay
#' details(res.dvar)
#'
#' # show coefficients estimated for each cell type
#' coefNames(res.dvar)
#'
#' # extract results using limma-style syntax
#' # combines all cell types together
#' # adj.P.Val gives study-wide FDR
#' topTable(res.dvar, coef = "group_idstim", number = 3)
#'
#' # Plot top hit to see differential variance
#' # Note that this is a toy example with only 4 samples
#' cellType <- 'CD4 T cells'     
#' gene <- 'DYNLRB1'
#' 
#' y <- res.proc[[cellType]]$E[gene,]
#' x <- colData(res.proc)$group_id
#' 
#' boxplot(y ~ x, 
#' 	xlab = "Stimulation status", 
#' 	ylab = "Gene expression",
#' 	main = paste(cellType, gene))
#
#' @seealso \code{variancePartition::diffVar()}, \code{missMethyl::diffVar()} 
#' @rdname diffVar-methods
#' @aliases diffVar,dreamletResult,dreamletResult-method
setMethod("diffVar", "dreamletResult",
	function( fit, method = c("AD", "SQ"),
		scale = c("leverage", "none"),
		BPPARAM = SerialParam(), ...){

	# run diffVar for each cell type
	fitList = lapply(fit, function(x){

		# test if any hatvalues are 1
		if( any(abs(x$hatvalues - 1.0) <= .Machine$double.eps) ){
			result = NULL
		}else{
			result = diffVar(x, method=method, scale = scale, BPPARAM=BPPARAM,... )
		}
		result
		})
	
	# keep only results that are not NULL
	keep = !sapply(fitList, is.null)
	fitList = fitList[keep]

	df_details = details(fit)
	i = match(names(fitList), df_details$assay)
	df_details = df_details[i,]

	# store results as dreamletResult
	new("dreamletResult", fitList, df_details = df_details)
})

# fit = readRDS("fit.RDS")
# fit2 = diffVar(fit)

# tab = topTable(fit, coef="c15xAD - c15xControl", number=Inf)
# tab2 = topTable(fit2, coef="c15xAD - c15xControl", number=Inf)

# df = merge(tab, tab2, by=c("assay", "ID"))
# df = as.data.frame(df)

# fig = ggplot(df, aes(t.x, t.y)) +
# 	geom_point() +
# 	theme_classic() +
# 	theme(aspect.ratio=1) +
# 	facet_wrap(~assay) +
# 	xlab("Differential expression") +
# 	ylab("Differential variance") +
# 	geom_smooth(color="red", method="lm")

# ggsave(fig, file="~/www/test.png")



# file = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/processAssays/MSSM_2023-09-12_17_04_processAssays_SubID_subclass.RDS"

# res.proc = readRDS( file )


# CT = 'Micro' 
# gene = 'IL15'

# CT = 'Oligo' 
# gene = 'NAALADL2'

# y = assay(res.proc, CT)$E[gene,]
# x = colData(res.proc)[,'AD',drop=FALSE]
# id = intersect(names(y), rownames(x))


# df2 = data.frame(x = as.factor(x[id,]), y=y[id])

# fig = ggplot(df2, aes(x, y)) +
# 	geom_violin() +
# 	geom_boxplot(width=.1) +
# 	theme_classic() +
# 	theme(aspect.ratio=1) +
# 	ggtitle(gene) +
# 	ylab(bquote(log[2]~Gene~expression)) +
# 	xlab("Alzheimer's disease status")

# ggsave(fig, file="~/www/test.pdf", width=5, height=5)

# tab[tab$ID == "NAALADL2",]

