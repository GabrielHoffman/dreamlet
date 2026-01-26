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
#' # Replace space with underscore, and remove "+" to avoid issue downstream
#' example_sce$cluster_id <- gsub(" ", "_", example_sce$cluster_id)
#' example_sce$cluster_id <- gsub("\\+", "", example_sce$cluster_id)
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
#' pb.stack <- stackAssays(pb)
#'
#' # voom-style normalization
#' # stackedAssay (i.e. cell type) can now be included as a covariate
#' res.proc <- processAssays(pb.stack, ~ group_id + stackedAssay)
#'
#' # Examine coding of covariates
#' # colData:
#' head(colData(res.proc))
#'
#' # Examine coding of covariates
#' # metadata:
#' head(metadata(res.proc))
#'
#' # Variance partitioning analysis
#' # Model contribution of Donor (id), stimulation status (group_id), and cell type (stackedAssay)
#' form <- ~ (1|id) + (1|group_id) + (1|stackedAssay)
#' vp <- fitVarPart(res.proc, form)
#'
#' # Summarize variance fractions across cell types
#' plotVarPart(sortCols(vp))
#'
#' # Interaction analysis allows group_id
#' # to have a different effect within each stackedAssay
#' form <- ~ (1|id) + (1|group_id) + (1|stackedAssay) + (1|group_id:stackedAssay)
#' vp2 <- fitVarPart(res.proc, form)
#'
#' plotVarPart(sortCols(vp2))
#' 
#' plotVarPart(sortCols(vp2))
#' 
#' # Differential expression analysis
#' # Testing differences between cell types
#' 
#' # In a real data you want to test the full model, 
#' # but this dataset is too small 
#' form <- ~ (1|id) + (1|group_id) + (1|stackedAssay) + group_id:stackedAssay + 0 
#'
#' # In this small dataset, just test simulation-by-celltype interaction term
#' # Here, test if the effect of stimulation (i.e. difference between stimulated and 
#' # controls) is different between B cells and monocytes
#' contrasts <- c(Diff = "(group_idstim:stackedAssayB_cells - group_idctrl:stackedAssayB_cells) - 
#' (group_idstim:stackedAssayCD14_Monocytes - group_idstim:stackedAssayCD14_Monocytes)")
#' form <- ~ group_id:stackedAssay + 0 
#' fit <- dreamlet( res.proc, form, contrasts = contrasts)
#'
#' # Top genes
#' topTable(fit, coef='Diff', number=3)
#' 
#' # Plot example
#' df <- extractData(res.proc, assay = "stacked", genes = c("ISG20"))
#' 
#' df <- df[df$stackedAssay %in% c("B_cells", "CD14_Monocytes"),]
#' 
#' ggplot(df, aes(group_id, ISG20)) +
#'   geom_boxplot() +
#'   theme_bw() +
#'   theme(aspect.ratio=1) +
#'   facet_wrap( ~ stackedAssay) 
#
#' @importFrom SummarizedExperiment assayNames<-
#' @importFrom rlang sym
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom dplyr rename inner_join
#' @export
stackAssays <- function(pb, assays = assayNames(pb)) {
  stopifnot(is(pb, "SingleCellExperiment"))

  # cbind assays
  counts <- lapply(assays, function(x) {
    y <- assay(pb, x)

    colnames(y) <- paste(x, colnames(y), sep = "_")
    y
  })
  counts <- do.call(cbind, counts)

  # rbind colData
  info <- lapply(assays, function(x) {
    info <- colData(pb)
    info$stackedAssay <- x
    rownames(info) <- paste(x, rownames(info), sep = "_")

    info
  })
  info <- do.call(rbind, info)

  # add metadata?

  # cell counts
  ids <- names(int_colData(pb)$n_cells)
  grd <- expand.grid(id = ids, assay = assays)
  rownames(grd) <- paste(grd$assay, grd$id, sep = "_")

  ncell.lst <- lapply(seq(nrow(grd)), function(i) {
    id <- as.character(grd$id[i])
    ct <- as.character(grd$assay[i])

    count <- int_colData(pb)$n_cells[[id]][ct]

    names(count) <- "stacked"
    count
  })
  names(ncell.lst) <- rownames(grd)

  # create SingleCellExperiment
  pb.stack <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = info
  )
  assayNames(pb.stack) <- "stacked"

  # include cell count data
  int_colData(pb.stack)$n_cells <- ncell.lst[colnames(pb.stack)]

  # metadata
  # agg_pars
  metadata(pb.stack)$agg_pars <- metadata(pb)$agg_pars

  # aggr_means
  # create structure, but don't populate values
  df <- grd[colnames(pb.stack), , drop = FALSE] %>%
    as_tibble()
  key <- metadata(pb.stack)$agg_pars$by
  df[[key[1]]] <- "stacked"
  df[[key[2]]] <- rownames(grd)

  # get metadata, joint with df
  by <- metadata(pb)$agg_pars$by
  df2 <- metadata(pb)$aggr_means %>%
    rename(
      assay = sym(by[1]),
      id = sym(by[2])
    )

  metadata(pb.stack)$aggr_means <- inner_join(df, df2, by = c("assay", "id"))

  pb.stack
}
