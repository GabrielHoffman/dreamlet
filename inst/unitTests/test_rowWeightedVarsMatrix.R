


test_rowWeightedVarsMatrix = function(){
		
	set.seed(1)
	n = 100
	p = 200
	X = matrix(rnorm(n*p), n, p)
	W = matrix(rgamma(n*p, 14, 1), n, p)

	# adapted from modi::weighted.var()0
	# placed here to avoid importing the whole package
	wv = function(x, w){
		w <- w * length(w)/sum(w)
    	sum(w * (x - weighted.mean(x, w))^2)/(sum(w) - 1)
	}

	res1 <- sapply(seq(n), function(i) wv(X[i,], W[i,]))

	res2 <- dreamlet:::rowWeightedVarsMatrix(X, W)

	checkEqualsNumeric(res1, res2)
}






