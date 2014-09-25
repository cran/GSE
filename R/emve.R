emve <- function(x, n.resample=100, maxits=5, n.sub.size)
{
	xcall <- match.call()

	## check dat
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")
	if( ncol(x) < 2 ) stop("Column dimension of 'x' must be at least 2.")

	
	## Check number of resampling
	if( n.resample < 20 ) stop("Number of resampling must be >= 20.")

	## drop all rows with missing values (!!) :
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	
	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")

	## Remove any rows with completely missing entries
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	if( length(not.ok) > 0) cat("Observations (rows): ", paste(not.ok, collapse=", "), 
		"\nare completely missing and will be dropped out from the estimation.\n")
	x_orig <- x
	x <- x[ ok,]
	x_nonmiss <- x_nonmiss[ ok,]
	
	## reorder the data based on missingness
	x_sort <- .sort.missing(x, x_nonmiss)

	## Check dimension 
	n <- nrow(x_sort$x); p <- ncol(x_sort$x)
	if(n <= p + 1)
		stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
	if(n < 2 * p)
		## p+1 < n < 2p
		warning("n < 2 * p, i.e., possibly too small sample size")

	## get initial estimate for all the EM calculation in EMVE subsampling
	mu0 <- colMeans(x_sort$x, na.rm=T)
	S0 <- diag(apply(x_sort$x, 2, var, na.rm=T))
	x_sort <- c(x_sort, .CovEM.setparam(p, mu0, S0))
	
	## Check subsample size
	if( missing(n.sub.size)) 
		n.sub.size <- floor( (p+1)/(1-mean(is.na(x))) )+1   #initial sub-sampling size
	if (n.sub.size <= p) 
		stop("n.sub.size <= p -- you can't be serious!")
			
	res <- with(x_sort, .emve.init(x, x_nonmiss, pu, n, p, theta, G.ind-1, length(theta),x.miss.group.match,
		miss.group.unique, miss.group.counts, miss.group.obs.col, miss.group.mis.col, 
		miss.group.p, miss.group.n, n.resample, maxits, n.sub.size))


	S.chol <- tryCatch( chol(res$S), error=function(e) NA)
	if( !is.matrix(S.chol) )  warning("Estimated covariance matrix is not positive definite. May consider increase the sample size of the data.")		

	pmd <- pmd.adj <- rep(NA, nrow(x_orig))
	pmd[ok] <- res$pmd[x_sort$id.ro]
	pmd.adj[ok] <- res$pmd.adj[x_sort$id.ro]
	pu <- rowSums( !is.na(x_orig))
	res <- new("emve",
		call = xcall,
		S = res$S,
		mu = res$mu,
		sc=res$mve.scale,
		cand.sc=res$cand.mve.scale,
		cand.mu=res$cand.mu,
		cand.S=res$cand.Sigma,
		estimator = "Extended Minimum Volume Ellipsoid", 
		x = x_orig,
		pmd = pmd,
		pmd.adj = pmd.adj,
		p = p, 
		pu = pu)
}
	
	
	
	
## Assume the input data matrix is sorted using .sort.missing
.emve.init <- function(x, x_nonmiss, pu, n, p, theta0, G, d, x.miss.group.match, miss.group.unique, miss.group.counts, 
	miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, n.resample=100, maxits=5, n.sub.size)
{
	##########################################################################################
	## SOME VARIABLE SETUP
	## Check subsample size
	if( missing(n.sub.size)) 
		n.sub.size <- floor( (p+1)/(1-mean(is.na(x))) )+1   #initial sub-sampling size
	if (n.sub.size <= p) 
		stop("n.sub.size <= p -- you can't be serious!")
		
	## Initial scales for later calculation
	scale0_init <- .scale.mve.init(pu)  # will be used to calculate EMVE scale
	cc <- scale0_init$cc
	ck <- scale0_init$ck
	colmed <- apply(x, 2, median, na.rm=T)	# coordinate-wise medians for each column, full data

	## Filled each entry by the overall median for that variable
	x_filled <- x
	a <- is.na(x)%*%diag(colmed)
	x_filled[is.na(x_filled)] <- 0
	x_filled <- x_filled + a
	##########################################################################################
	
	res <- .emve.Rcpp(x_filled, x_nonmiss, pu, n, p, theta0, G, d, x.miss.group.match, miss.group.unique, miss.group.counts,
			miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, n.resample, n.sub.size, 0.00001, cc, ck, maxits)

	## only the best 10 candidates
	res <- list(
		S = res$S0[,,1],
		mu = res$mu0[1,],
		mve.scale=res$ss0[1],
		cand.mve.scale=res$ss0[1:10],
		cand.mu=res$mu0[1:10,],
		cand.Sigma=res$S0[,,1:10])
		
	x.tmp <- sweep(x, 2, res$mu, "-")	
	pmd <- .partial.mahalanobis.Rcpp(x.tmp, res$S, miss.group.unique, miss.group.counts)
	
	## compute adjusted pmd
	pmd.adj <- qchisq( pchisq( pmd, df=pu, log.p=T, lower.tail=F), df=p, log.p=T, lower.tail=F) 
	pmd.adj[ which( pu == p)] <- pmd[ which(pu==p) ]

	res$pmd <- pmd
	res$pmd.adj <- pmd.adj
	res
}

	
	
	


