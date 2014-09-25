.gy.filt <- function(x, alpha, it=TRUE){
	x <- as.matrix(x)
	xs <- scale(x, apply(x, 2, median, na.rm=TRUE), apply(x, 2, mad, na.rm=TRUE))
	xs2 <- xs^2
	xs2.na <- switch( it + 1, 
		apply(xs2, 2, .gy.filt.uni, alpha=alpha),
		apply(xs2, 2, .gy.filt.uni.it, alpha=alpha))
	x.na <- x
	x.na[ which(is.na(xs2.na)) ] <- NA
	if( ncol(x.na) == 1) x.na <- c(x.na)
	return( x.na )
}


.gy.filt.uni <- function(v, alpha){
	n <- length(v)
	v.order <- order(v)	
	v <- sort(v)
	i0 <- which(v < qchisq( alpha, 1 ))
	n0 <- 0
	if( length(i0) > 0){
		i0 <- rev(i0)[1]
		dn <- max( pmax( pchisq( v[i0:n], 1) - (i0:n - 1)/n, 0)) 
		n0 <- round(dn*n)
	} 
	v <- v[ order(v.order) ] 
	v.na <- v
	if(n0 > 0) v.na[ v.order[ (n - n0 + 1):n] ] <- NA
	return( v.na )
}

.gy.filt.uni.it <- function(v, alpha, miter=30){
	converge <- 0	
	iter <- 0
	n <- length(v)
	id <- 1:n
	v.old <- v
	while( converge == 0 & iter < miter ){
		iter <- iter + 1
		v <- .gy.filt.uni( v, alpha )
		id <- id[ !is.na(v) ] 
		if( !any(is.na(v)) ) converge <- 1
		v <- na.omit(v)
	}
	v.out <- rep(NA, n)
	v.out[id] <- v
	return( v.out)
}


TSGS <- function(x, alpha=0.99, it=TRUE, ...){
	xcall <- match.call()
	
	## check dat
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame.")
	if(any(is.na(x))) warning("Data matrix contains missing values.")
	## June 12, 2012
	## Only allow up to p=50
	p <- ncol(x)
	if( p >200 | p < 2 ) stop("Column dimension of 'x' must be in between 2 and 200.")
	
	xf <- .gy.filt(x, alpha=alpha, it=it)
	res <- GSE(xf, ...)
	res <- new("TSGS",
		call = xcall,
		S = res@S,
		mu = res@mu,
		xf = xf,
		sc = res@sc,
		mu0 = res@mu0,
		S0 = res@S0, 
		iter = res@iter,
		eps = res@eps,
		estimator = "2SGS", 
		x = x,
		ximp = res@ximp,
		weights = res@weights,
		weightsp = res@weightsp,
		pmd = res@pmd,
		pmd.adj = res@pmd.adj,
		p = res@p,
		pu = res@pu)
	res
}
