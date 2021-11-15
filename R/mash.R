# # Gabriel Hoffman
# # Nov 15, 2021


# #' Convert results table to matrix
# #' 
# #' Convert results table to matrix
# #' 
# #' @param tab results table from \code{topTable()}
# #' @param col which column to extract
# #' 
# #' @importFrom Matrix sparseMatrix 
# toMatrix = function(tab, col){
  
#   # row and column names
#   rn = unique(tab$ID)
#   cn = unique(tab$assay)

#   i = match(tab$ID, rn)
#   j = match(tab$assay, cn)

#   M = sparseMatrix(i,j, x=tab[[col]], 
#     dims=c(length(rn), length(cn)),
#     dimnames = list(rn, cn))

#   data = as.matrix(M)
#   data[data == 0] = NA 
#   data
# }

# # fit = res.dl
# # coef = 'StimStatusstim'

# #' @importFrom mashr mash_set_data cov_canonical mash_estimate_corr_em
# run_mash = function( fit, coef){

# 	# get results for each gene and cell type
# 	tab = topTable(fit, coef=coef, Inf)

# 	# compute standard error from t-stat and logFC
# 	tab$se = tab$logFC / tab$t

# 	# convert to matricies
# 	B = toMatrix(tab, "logFC")
# 	S = toMatrix(tab, "se")

# 	# run mashr on these matricies
# 	#-----------------------------

# 	# set up
# 	# NA's are replaced with beta = 0 with se = 1e6 
# 	data = mash_set_data(B, S)

# 	# estimate some parameters
# 	U.c = cov_canonical(data)

# 	# Estimate correlation structure 
# 	V.em = mash_estimate_corr_em(data, U.c, details = TRUE)

# 	# get model fit
# 	V.em$mash.model
# }


# # see posterior mean for logFC
# # head(get_pm(m.Vem))

# # # how many gene-by-celltype are significant
# # # i.e.  if a gene is significant in 2 celltypes, it is counted twice
# # table(get_lfsr(m.Vem) < 0.05)

# # # how many genes are significant in at least one cell type
# # table( apply(get_lfsr(m.Vem), 1, min) < 0.05)

# # # how many genes are significant in each cell type
# # apply(get_lfsr(m.Vem), 2, function(x) sum(x < 0.05))
# # ```

# # # examine top set of genes
# # ```{r mashr_summary2}
# # # which genes are significant in at least 1 cell type
# # sort(names(get_significant_results(m.Vem)))

# # # Lets examine APOE.
# # # There is a lot of variation in the raw logFC
# # B["APOE",]

# # # posterior mean after borrowing across cell type and genes
# # # There might be too much borrowing!
# # get_pm(m.Vem)["APOE",]




