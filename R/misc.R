.sort.missing <- function(x_nonmiss){
	## Sort the observations based on their missingness for faster computation (in partial mahalanobis distance)
	## Added Mar 17, 2012
	miss.group <- apply( x_nonmiss, 1, paste, collapse='')
	id.order <- order(miss.group)
	miss.group.unique <- do.call(rbind, strsplit(sort(unique(miss.group)),""))
	miss.group.unique <- t(apply(miss.group.unique, 1, as.numeric))
	miss.group.counts <- as.numeric(table(miss.group))

	return( list(id.order=id.order, miss.group.unique=miss.group.unique, 
		miss.group.counts=miss.group.counts) )
}


## Compute partial mahalanobis distance 
## this is for the case when the data set is not sorted 
partial.mahalanobis <- function(x, mu, Sigma){
	xcall <- match.call()

	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")

	## drop all rows with missing values (!!) :
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	x.orig <- x
	x <- x[ ok,]
	x_nonmiss <- x_nonmiss[ ok,]
	n <- nrow(x); p <- ncol(x) 

	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")
	
	if( any( x_nonmiss == 0 ) ){
		x.sort.miss <- .sort.missing(x_nonmiss)
		id.order <- x.sort.miss$id.order
		x <- x[id.order,]
		x_nonmiss <- x_nonmiss[id.order,]
		pp <- pp[id.order]	
		miss.group.unique <- x.sort.miss$miss.group.unique
		miss.group.counts <- x.sort.miss$miss.group.counts
	} else{
		miss.group.unique <- matrix( 1, nrow=1, ncol=p )
		miss.group.counts <- n
		id.order <- 1:n
	}
	x <- sweep(x, 2, c(mu), "-")
	pmd.tmp <- .partial.mahalanobis.Rcpp( x, Sigma, miss.group.unique, miss.group.counts )
	pmd.tmp <- pmd.tmp[ order(id.order) ]
	
	pmd <- rep(NA, nrow(x.orig))
	pmd[ok] <- pmd.tmp
	
	pu <- rowSums( !is.na(x.orig))
	pmd.adj <- qchisq( pchisq( pmd, df=pu, log.p=T, lower.tail=F), df=p, log.p=T, lower.tail=F) 
	pmd.adj[ which( pu == p)] <- pmd[ which(pu==p) ]	
	
	new("CovRobMiss", 
		mu = mu,
		S = Sigma,
		call = xcall,
		estimator = "unknown",
		x = x.orig, 
		pmd = pmd,
		pmd.adj = pmd.adj,
		p=p,
		pu = pu)	
}

