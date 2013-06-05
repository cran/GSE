###########################################################
## Generalized S-Estimator
###########################################################
GSE <- function(x, tol=1e-5, maxiter=500, init="emve", tol.scale=1e-4, miter.scale=30, print.step=1, mu0, S0, ...)
{
	## argument checks
	## check choices of psi
	init <- match.arg(tolower(init), choices=c("emve","sign","qc","huber") )
	
	if( !is.numeric(print.step) | print.step < 0 | print.step > 2) stop("argument 'print.step' must be: 0, 1, 2.")
	
	## check dat
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")

	## drop all rows with missing values (!!) :
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	if( length(not.ok) > 0) cat("Observations (rows): ", paste(not.ok, collapse=", "), 
		"\nare completely missing and will be dropped out from the estimation.\n")
	x.orig <- x
	x <- x[ ok,]
	x_nonmiss <- x_nonmiss[ ok,]

	## reorder the data based on missingness
	x.sort.miss <- .sort.missing(x_nonmiss)
	miss.group.unique <- x.sort.miss$miss.group.unique
	miss.group.counts <- x.sort.miss$miss.group.counts
	id.order <- x.sort.miss$id.order
	
	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")

	## dimension
	n <- nrow(x); p <- ncol(x)
	if(n <= p + 1)
		stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
	if(n < 2 * p)
		## p+1 < n < 2p
		warning("n < 2 * p, i.e., possibly too small sample size")

	xcall <- match.call()

	if( xor(missing(mu0), missing(S0)) ) warning("Both 'mu0' and 'S0' must be provided. Default 'init' is used...")
	if( missing(mu0) || missing(S0) ){
		init.res <- switch( init,
			emve=.emve.init(x, miss.group.unique, miss.group.counts, id.order, ...),
			qc ={res <- HuberPairwise(x, psi="sign", computePmd = FALSE); list(mu=res@mu, S=res@S) },
			sign ={res <- HuberPairwise(x, psi="sign", computePmd = FALSE); list(mu=res@mu, S=res@S)},			
			huber = {res <- HuberPairwise(x, psi="huber", computePmd = FALSE, ...); list(mu=res@mu, S=res@S)} )
		S0 <- init.res$S
		mu0 <- init.res$mu		
	} 
	
	## initiate GSE computation
	res <- .GSE.init(x, miss.group.unique, miss.group.counts, id.order, mu0, S0, tol, maxiter, tol.scale, miter.scale, print.step)

	## compute pmd
	pmd <- rep(NA, nrow(x.orig))
	pmd[ok] <- res$pmd

	## compute adjusted pmd
	pmd.adj <- qchisq( pchisq( pmd, df=pp, log.p=T, lower.tail=F), df=p, log.p=T, lower.tail=F) 
	pmd.adj[ which( pp == p)] <- pmd[ which(pp==p) ]
	
	res <- new("GSE",
		call = xcall,
		S = res$S,
		mu = res$mu,
		sc = res$stilde0,
		iter = res$iter,
		eps = res$ep,
		estimator = "Generalized S-Estimator", 
		x = x.orig,
		pmd = pmd,
		pmd.adj = pmd.adj,
		p = p,
		pu = pp)
	res
}


.GSE.init <- function(x, miss.group.unique, miss.group.counts, id.order, mu0, S0, tol, maxiter, tol.scale, miter.scale, print.step)
{
	##########################################################################################
	## basic variables initialization
	n <- nrow(x)
	p <- ncol(x)
	x <- x[id.order,]
	tuning.const.group <- .rho.tune(apply(miss.group.unique,1,sum))

	## June 12, 2012
	## Only allow up to p=50
	if( p>200) stop("Current version only allows p <= 200.")
	
	##########################################################################################
	## Start computing
	res <- .GSE.Rcpp(x, matrix(mu0,1,p), S0, tol, maxiter, tol.scale, miter.scale, 
		miss.group.unique, miss.group.counts, tuning.const.group, print.step)
	
	##########################################################################################
	## Include additional output other than mu and S output by .GSE.fixpt.init:
	## input data 'x', nonmissing variables 'pu', partial MD 'pmd', 
	## note: need to reorder x, pmd, pu
	id.reorder <- order(id.order)
	
	## Compute pmd
	x.tmp <- sweep(x, 2, c(res$mu), "-")
	res$pmd <- .partial.mahalanobis.Rcpp(x.tmp, res$S, miss.group.unique, miss.group.counts)
	res$pmd <- res$pmd[ id.reorder]

	##########################################################################################
	## Output
	res
}







