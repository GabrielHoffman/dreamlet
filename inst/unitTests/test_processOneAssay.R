

test_processOneAssay = function(){

    data(varPartDEdata)
    
	set.seed(1)
	mix = sample.int(ncol(countMatrix), ncol(countMatrix))
	dge = DGEList(counts = countMatrix[,mix])
	dge = calcNormFactors(dge)
	rownames(metadata) = rownames(metadata)[mix]

	a = with(dge$samples, lib.size * norm.factors)
	b = getNormLibSizes(dge)
	checkEqualsNumeric(a,b)

	form <- ~ Disease 
	dsgn = model.matrix(form, metadata)

	# processOneAssay
	#------------------

	w = 1:ncol(dge)
	fit1r = voomLmFit2( dge, dsgn, prior.weights=w / mean(w))
	fit1 = eBayes(fit1r)

	wMat = matrix(w / mean(w), nrow = nrow(dge), ncol = length(w), byrow=T)
	rownames(wMat) = rownames(dge)
	colnames(wMat) = colnames(dge)

	res = dreamlet:::processOneAssay(dge, form, metadata, 
			weights = wMat, 
			n.cells = w,
			min.cells = 0, 
			min.count = 0, 
			min.samples = 0, 
            prior.count = 0.5,
			min.prop = 0, 
			span = 0.5)

	checkEquals(fit1$EList$E, res$E, tolerance = 1e-4)
	checkEquals(as.numeric(fit1$EList$weights), as.numeric(res$weights), tolerance = 3e-4)
}





# Fig bug reported here
# https://support.bioconductor.org/p/9154670/


voomLmFit2 = function (counts, design = NULL, block = NULL, prior.weights = NULL, 
    sample.weights = FALSE, var.design = NULL, var.group = NULL, 
    lib.size = NULL, normalize.method = "none", span = 0.5, plot = FALSE, 
    save.plot = FALSE, keep.EList = TRUE) 
{
    Block <- !is.null(block)
    PriorWeights <- !is.null(prior.weights)
    SampleWeights <- sample.weights || !is.null(var.design) || 
        !is.null(var.group)
    if (PriorWeights && SampleWeights) 
        stop("Can't specify prior.weights and estimate sample weights")
    out <- list()
    if (is(counts, "SummarizedExperiment")) 
        counts <- SE2DGEList(counts)
    if (is(counts, "DGEList")) {
        out$genes <- counts$genes
        out$targets <- counts$samples
        if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
            0) 
            design <- model.matrix(~group, data = counts$samples)
        if (is.null(lib.size)) 
            lib.size <- getNormLibSizes(counts)
        counts <- counts$counts
    }
    else {
        if (is(counts, "eSet")) {
            if (!requireNamespace("Biobase", quietly = TRUE)) 
                stop("Biobase package required but is not installed (or can't be loaded)")
            if (length(Biobase::fData(counts))) 
                out$genes <- Biobase::fData(counts)
            if (length(Biobase::pData(counts))) 
                out$targets <- Biobase::pData(counts)
            counts <- get("counts", Biobase::assayData(counts))
        }
        else {
            counts <- as.matrix(counts)
        }
    }
    n <- nrow(counts)
    if (n < 2L) 
        stop("Need at least two genes to fit a mean-variance trend")
    m <- min(counts)
    if (is.na(m)) 
        stop("NA counts not allowed")
    if (m < 0) 
        stop("Negative counts not allowed")
    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }
    if (is.null(lib.size)) 
        lib.size <- colSums(counts)
    if (!is.null(prior.weights)) 
        prior.weights <- asMatrixWeights(prior.weights, dim(counts))
    # added GEH
    attr(prior.weights, "arrayweights") = NULL
    y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
    y <- normalizeBetweenArrays(y, method = normalize.method)
    fit <- lmFit(y, design, weights = prior.weights)
    if (is.null(fit$qr)) 
        h <- hat(design, intercept = FALSE)
    else h <- hat(fit$qr)
    MinGroupSize <- 1/max(h)
    eps <- 1e-04
    RowHasZero <- which(rowSums(counts < eps) > (max(2, MinGroupSize) - 
        eps))
    AnyZeroRows <- as.logical(length(RowHasZero))
    if (AnyZeroRows) {
        countsZero <- counts[RowHasZero, , drop = FALSE]
        PoissonFit <- glmFit(countsZero, design = design, lib.size = lib.size, 
            dispersion = 0, prior.count = 0)
        IsZero <- (PoissonFit$fitted.values < eps & countsZero < 
            eps)
        RowHasExactZero <- which(rowSums(IsZero) > eps)
        if (length(RowHasExactZero)) {
            RowHasZero <- RowHasZero[RowHasExactZero]
            IsZero <- IsZero[RowHasExactZero, , drop = FALSE]
            yNAshort <- y[RowHasZero, , drop = FALSE]
            yNAshort[IsZero] <- NA
            fitNA <- suppressWarnings(lmFit(yNAshort, design, 
                weights = prior.weights[RowHasZero, , drop = FALSE]))
            fit$df.residual[RowHasZero] <- fitNA$df.residual
            fit$sigma[RowHasZero] <- fitNA$sigma
            if (Block || SampleWeights) {
                yNAfull <- y
                yNAfull[RowHasZero, ] <- yNAshort
            }
        }
        else {
            AnyZeroRows <- FALSE
        }
    }
    HasRep <- (fit$df.residual > 0L)
    NWithReps <- sum(HasRep)
    if (NWithReps < 2L) {
        if (NWithReps == 0L) 
            warning("The experimental design has no replication. Setting weights to 1.")
        if (NWithReps == 1L) 
            warning("Only one gene with any replication. Setting weights to 1.")
        fit$genes <- out$genes
        return(fit)
    }
    Amean <- Amean2 <- rowMeans(y)
    if (AnyZeroRows) 
        Amean2[RowHasZero] <- rowMeans(yNAshort, na.rm = TRUE)
    sx <- Amean2[HasRep] + mean(log2(lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma[HasRep])
    if (AnyZeroRows) 
        l <- weightedLowess(sx, sy, span = span, weights = fit$df.residual[HasRep], 
            output.style = "lowess")
    else l <- lowess(sx, sy, f = span)
    if (plot) {
        plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
            pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lty <- ifelse(Block || SampleWeights, 2, 1)
        lines(l, col = "red", lty = lty)
    }
    f <- approxfun(l, rule = 2, ties = list("ordered", mean))
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coefficients[, j, drop = FALSE] %*% 
            t(fit$design[, j, drop = FALSE])
    }
    else {
        fitted.values <- fit$coefficients %*% t(fit$design)
    }
    fitted.cpm <- 2^fitted.values
    fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    if (PriorWeights) 
        weights <- w * prior.weights
    else weights <- w

    # weights = weights / rowMeans(weights)
    
    if (SampleWeights) {
        if (AnyZeroRows) {
            sw <- arrayWeights(yNAfull, design, weights = weights, 
                var.design = var.design, var.group = var.group)
        }
        else {
            sw <- arrayWeights(y, design, weights = weights, 
                var.design = var.design, var.group = var.group)
        }
        message("First sample weights (min/max) ", paste(format(range(sw)), 
            collapse = "/"))
        if (Block) 
            weights <- t(sw * t(weights))
    }
    if (Block) {
        if (AnyZeroRows) {
            dc <- suppressWarnings(duplicateCorrelation(yNAfull, 
                design, block = block, weights = weights))
        }
        else {
            dc <- suppressWarnings(duplicateCorrelation(y, design, 
                block = block, weights = weights))
        }
        correlation <- dc$consensus.correlation
        if (is.na(correlation)) 
            correlation <- 0
        message("First intra-block correlation  ", format(correlation))
    }
    else {
        correlation <- NULL
    }
    if (Block || SampleWeights) {
        if (SampleWeights) 
            weights <- asMatrixWeights(sw, dim(y))
        else weights <- prior.weights
        fit <- lmFit(y, design, block = block, correlation = correlation, 
            weights = weights)
        if (AnyZeroRows) {
            fitNA <- suppressWarnings(lmFit(yNAshort, design, 
                block = block, correlation = correlation, weights = weights[RowHasZero, 
                  , drop = FALSE]))
            fit$df.residual[RowHasZero] <- fitNA$df.residual
            fit$sigma[RowHasZero] <- fitNA$sigma
        }
        sy <- sqrt(fit$sigma[HasRep])
        if (AnyZeroRows) 
            l <- weightedLowess(sx, sy, span = span, weights = fit$df.residual[HasRep], 
                output.style = "lowess")
        else l <- lowess(sx, sy, f = span)
        if (plot) {
            lines(l, col = "red")
            legend("topright", lty = c(2, 1), col = "red", legend = c("First", 
                "Final"))
        }
        f <- approxfun(l, rule = 2, ties = list("ordered", mean))
        if (fit$rank < ncol(design)) {
            j <- fit$pivot[1:fit$rank]
            fitted.values <- fit$coefficients[, j, drop = FALSE] %*% 
                t(fit$design[, j, drop = FALSE])
        }
        else {
            fitted.values <- fit$coefficients %*% t(fit$design)
        }
        fitted.cpm <- 2^fitted.values
        fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 
            1))
        fitted.logcount <- log2(fitted.count)
        w <- 1/f(fitted.logcount)^4
        dim(w) <- dim(fitted.logcount)
        if (PriorWeights) 
            weights <- w * prior.weights
        else weights <- w
        if (SampleWeights) {
            if (AnyZeroRows) {
                sw <- arrayWeights(yNAfull, design, weights = weights, 
                  var.design = var.design, var.group = var.group)
            }
            else {
                sw <- arrayWeights(y, design, weights = weights, 
                  var.design = var.design, var.group = var.group)
            }
            message("Final sample weights (min/max) ", paste(format(range(sw)), 
                collapse = "/"))
            weights <- t(sw * t(weights))
        }
        if (Block) {
            if (AnyZeroRows) {
                dc <- suppressWarnings(duplicateCorrelation(yNAfull, 
                  design, block = block, weights = weights))
            }
            else {
                dc <- suppressWarnings(duplicateCorrelation(y, 
                  design, block = block, weights = weights))
            }
            correlation <- dc$consensus.correlation
            if (is.na(correlation)) 
                correlation <- 0
            message("Final intra-block correlation  ", format(correlation))
        }
    }
    fit <- lmFit(y, design, block = block, correlation = correlation, 
        weights = weights)
    if (is.null(fit$Amean)) 
        fit$Amean <- Amean
    if (AnyZeroRows) {
        fitNA <- suppressWarnings(lmFit(yNAshort, design, block = block, 
            correlation = correlation, weights = weights[RowHasZero, 
                , drop = FALSE]))
        fit$df.residual[RowHasZero] <- fitNA$df.residual
        fit$sigma[RowHasZero] <- fitNA$sigma
    }
    fit$genes <- out$genes
    fit$targets <- out$targets
    if (is.null(fit$targets)) {
        fit$targets <- data.frame(lib.size = lib.size)
        row.names(fit$targets) <- colnames(y)
    }
    if (SampleWeights) 
        fit$targets$sample.weight <- sw
    if (save.plot) {
        fit$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )", 
            ylab = "Sqrt( standard deviation )")
        fit$voom.line <- l
    }
    if (keep.EList) {
        fit$EList <- new("EList", list(E = y, weights = weights, 
            genes = out$genes))
    }
    fit
}












