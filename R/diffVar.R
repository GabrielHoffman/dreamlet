

# # @export
# # @rdname diffVar-method
# # @aliases diffVar,dreamletResult-method
# setMethod("diffVar", "dreamletResult",
# 	function( fit, method = c("AD", "SQ"),
# 		BPPARAM = SerialParam(), ...){

# 	# run diffVar for each cell type
# 	fitList = lapply(fit, function(x){

# 		# test if any hatvalues are 1
# 		if( any(abs(x$hatvalues - 1.0) <= .Machine$double.eps) ){
# 			result = NULL
# 		}else{
# 			result = diffVar(x, method=method, BPPARAM=BPPARAM,... )
# 		}
# 		result
# 		})
	
# 	# keep only results that are not NULL
# 	keep = !sapply(fitList, is.null)
# 	fitList = fitList[keep]

# 	df_details = details(fit)
# 	i = match(names(fitList), df_details$assay)
# 	df_details = df_details[i,]

# 	# store results as dreamletResult
# 	new("dreamletResult", fitList, df_details = df_details)
# })

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

