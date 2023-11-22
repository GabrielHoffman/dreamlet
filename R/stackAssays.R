# Gabriel Hoffman
# Nov 16, 2023

#' Stack assays from pseudobulk 
#'
#' Stack assays from pseudobulk to perform analysis across cell types
#' 
#' @param pb pseudobulk \code{SingleCellExperiment} from \code{aggregateToPseudoBulk()}
#' @param assays array of assay names to include in analysis. Defaults to \code{assayNames(pb)}
#' 
#' @return pseudobulk \code{SingleCellExperiment} cbind'ing expression values and rbind'ing colData.  The column \code{stackedAssay} in \code{colData()} stores the assay information of the stacked data.
#'
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
#' # Stack assays for joint analysis
#' pb.stack = stackAssays(pb)
#' 
#' # voom-style normalization
#' # assay (i.e. cell type) can now be included as a covariate
#' res.proc <- processAssays(pb.stack, ~ group_id + stackedAssay)
#'
#' # variance partitioning analysis
#' vp <- fitVarPart(res.proc, ~ group_id + stackedAssay)
#'     
#' # Summarize variance fractions across cell types
#' plotVarPart( sortCols( vp ) )
#'
#' # Interaction analysis allows group_id
#' # to have a different effect within each stacedAssay
#' vp2 <- fitVarPart(res.proc, ~ group_id*stackedAssay)
#'     
#' plotVarPart( sortCols( vp2 ) )
#'
#' # Interaction model using random effects
#' form <- ~ (1|group_id) + (1|stackedAssay) + (1|group_id:stackedAssay)
#' @importFrom SummarizedExperiment assayNames<-
#' @export
stackAssays <- function(pb, assays = assayNames(pb)){

	stopifnot( is(pb, "SingleCellExperiment") )

	# cbind assays
	counts <- lapply(assays, function(x){
		y <- assay(pb, x)

		colnames(y) <- paste(x, colnames(y), sep="_")
		y 
		})
	counts <- do.call(cbind, counts)

	# rbind colData
	info <- lapply(assays, function(x){
		info <- colData(pb)
		info$stackedAssay <- x
		rownames(info) <- paste(x, rownames(info), sep='_') 

		info
		})
	info <- do.call(rbind, info)

	# add metadata?

	# cell counts
	ids <- names(int_colData(pb)$n_cells)
	grd <- expand.grid( assay = assays, id = ids)
	rownames(grd) <- paste(grd$assay, grd$id, sep='_')

	ncell.lst <- lapply(seq(nrow(grd)), function(i){

		id <- as.character(grd$id[i])
		ct <- as.character(grd$assay[i])

		count <- int_colData(pb)$n_cells[[id]][ct]

		names(count) <- "stacked"
		count
		})
	names(ncell.lst) <- rownames(grd)

	# create SingleCellExperiment
	pb.stack <- SingleCellExperiment(assays = list(counts = counts),
						colData = info)
	assayNames(pb.stack) <- 'stacked'

	# include cell count data
	int_colData(pb.stack)$n_cells <- ncell.lst[colnames(pb.stack)]

	# metadata
	# agg_pars
	metadata(pb.stack)$agg_pars <- metadata(pb)$agg_pars

	# aggr_means
	# create structure, but don't populate values
	df <- grd[colnames(pb.stack),,drop=FALSE]
	key <- metadata(pb.stack)$agg_pars$by
	df[[key[1]]] <- "stacked"
	df[[key[2]]] <- rownames(grd)

	metadata(pb.stack)$aggr_means <- tibble(df[,key,drop=FALSE])

	pb.stack
}		








