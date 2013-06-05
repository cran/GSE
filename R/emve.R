emve <- function(x, n.resample=100, maxits=5, n.sub.size)
{
	## check dat
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")

	## Check number of resampling
	if( n.resample < 20 ) stop("Number of resampling must be >= 20.")

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

	## Check dimension 
	n <- nrow(x)
	p <- ncol(x)
	if(n <= p + 1)
		stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
	if(n < 2 * p)
		## p+1 < n < 2p
		warning("n < 2 * p, i.e., possibly too small sample size")

	## Check subsample size
	if( missing(n.sub.size)) 
		n.sub.size <- floor( (p+1)/(1-mean(is.na(x))) )+1   #initial sub-sampling size
	if (n.sub.size <= p) 
		stop("n.sub.size <= p -- you can't be serious!")
		
	xcall <- match.call()
	
	res <- .emve.init(x, miss.group.unique, miss.group.counts, id.order, n.resample, maxits, n.sub.size)
	
	S.chol <- tryCatch( chol(res$S), error=function(e) NA)
	if( !is.matrix(S.chol) )  warning("Estimated covariance matrix is not positive definite. May consider increase the sample size of the data.")		

	pmd <- pmd.adj <- rep(NA, nrow(x.orig))
	pmd[ok] <- res$pmd
	pmd.adj[ok] <- res$pmd.adj
	pu <- rowSums( !is.na(x.orig))
	res <- new("emve",
		call = xcall,
		S = res$S,
		mu = res$mu,
		sc=res$mve.scale,
		cand.sc=res$cand.mve.scale,
		cand.mu=res$cand.mu,
		cand.S=res$cand.Sigma,
		estimator = "Extended Minimum Volume Ellipsoid", 
		x = x.orig,
		pmd = pmd,
		pmd.adj = pmd.adj,
		p = p, 
		pu = pu)
}
	
	
.emve.init <- function(x, miss.group.unique, miss.group.counts, id.order, n.resample=100, maxits=5, n.sub.size)
{
	##########################################################################################
	## SOME VARIABLE SETUP
	n.candidates <- 20				# number of intermediate candidates
	nf <- 10						# number of final candidates
	p <- ncol(x)					# number of variable
	n <- nrow(x)					# sample size	
	
	##########################################################################################
	## Check subsample size
	if( missing(n.sub.size)) 
		n.sub.size <- floor( (p+1)/(1-mean(is.na(x))) )+1   #initial sub-sampling size
	if (n.sub.size <= p) 
		stop("n.sub.size <= p -- you can't be serious!")
		
	##########################################################################################
	## SOME DATA FORMATTING
	## Nonmissing entry matrix for x
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)

	## Initial estimate
	scale0_init <- .scale.mve.init(pp)  # will be used to calculate EMVE scale
	med <- apply(x, 2, median, na.rm=T)					# coordinate-wise medians for each column, full data

	## Store results
	C0 <- array(NA, c(p, p, n.candidates))
	c0 <- array(NA, c(n.candidates, p))
	ss0 <- rep(NA, n.candidates)

	## Filled each entry by the overall median for that variable
	x.filled <- x
	a <- is.na(x)%*%diag(med)
	x.filled[is.na(x.filled)] <- 0
	x.filled <- x.filled + a
	
	## Sort the observations based on their missingness for faster computation (in partial mahalanobis distance)
	## Added Mar 17, 2012
	x <- x[id.order,]
	x.filled <- x.filled[id.order,]
	x_nonmiss <- x_nonmiss[id.order,]
	pp <- pp[id.order]
	##########################################################################################
	
	##########################################################################################
	## INITIALIZATION, i.e. search for subsample with good conditioned covariance matrix 
	curMinRcond <- 0.00001
	x.res.init <- .emve.resamp(x.filled, x_nonmiss, n.resample, n.sub.size, curMinRcond)
	x.res.ind <- x.res.init$subsample_index_list
	curMinRcond <- x.res.init$curMinRcond
    rm(x.res.init)
	## END INITIALIZATION 	 
	########################################################################################## 
	
	##########################################################################################
	## START
	for(i in 1:n.resample)
	{
		numCandFilled <- length(na.omit(ss0))

		## Compute coordinate-wise median for the given subsample
		x.sub <- x[x.res.ind[i,],]
		x.sub.filled <- x.filled[x.res.ind[i,],]
		mu <- apply(x.sub, 2, median, na.rm=T)
		## Compute sample covariance of the completed subsample 
		## At this stage, all sample covariances should be in good conditioned
		Sigma <- cov(x.sub.filled)

		## Compute EMVE scale
		## First find a based on the sample covariance such that the constraint is satisfied
		a <- .emve.scale.missing.Rcpp(Sigma, miss.group.unique, miss.group.counts )$a
		Sigma <- a*Sigma

		x_mu_diff <- sweep(x.filled, 2, mu, "-")
		partial <- .partial.mahalanobis.Rcpp(x_mu_diff, Sigma, miss.group.unique, miss.group.counts)
		scale0 <- .fast.scale.mve(partial, scale0_init$cc, scale0_init$ck)
		
		if( numCandFilled < n.candidates ){
			ss0[ numCandFilled+1 ] <- scale0
			c0[ numCandFilled+1, ] <- mu
			C0[ , , numCandFilled+1] <- Sigma
		} else{
			curMaxScaleInd <- which.max( ss0 )
			if( scale0 < max(ss0, na.rm=T) ){
				ss0[ curMaxScaleInd ] <- scale0
				c0[ curMaxScaleInd, ] <- mu
				C0[ , , curMaxScaleInd] <- Sigma
			}				
		}

		## Concentration step: choose 50% of the pts with smallest 
		## pi_i or corrected distance and compute Gaussian MLE for that
		partial.corrected <- pchisq(partial/scale0,pp)
		res <- rank(partial.corrected)
		mm <- (1:n)[res<=floor(n/2)]
		x.half <- x[mm,]
		gaussian.mle <- .cov.EM.rough(x.half, maxits)
		##gaussian.mle <- cov.EM(x.half)
		mu <- gaussian.mle$mean
		Sigma <- gaussian.mle$cov
		rm(gaussian.mle)

		## Again compute the EMVE scale based on the Gaussian MLE
		if( is.matrix(Sigma) ){
			if( rcond(Sigma) >= curMinRcond ){
				a <- .emve.scale.missing.Rcpp( Sigma, miss.group.unique, miss.group.counts )$a
	      		Sigma <- a*Sigma

				x_mu_diff <- sweep(x.filled, 2, mu, "-")
				partial <- .partial.mahalanobis.Rcpp(x_mu_diff, Sigma, miss.group.unique, miss.group.counts)
				scale0 <- .fast.scale.mve(partial, scale0_init$cc, scale0_init$ck)
				
				numCandFilled <- length(na.omit(ss0))
				if( numCandFilled < n.candidates ){
					ss0[ numCandFilled+1 ] <- scale0
					c0[ numCandFilled+1, ] <- mu
					C0[ , , numCandFilled+1] <- Sigma
				} else{
					curMaxScaleInd <- which.max( ss0 )
					if( scale0 < max(ss0, na.rm=T) ){
						ss0[ curMaxScaleInd ] <- scale0
						c0[ curMaxScaleInd, ] <- mu
						C0[ , , curMaxScaleInd] <- Sigma
					}				
				}
			}
		}
	}

	## Compute Gaussian MLE based on the entire sample
	gaussian.mle <- .cov.EM.rough(x, maxits*2)
	##gaussian.mle <- cov.EM(x)
	mu <- gaussian.mle$mean
	Sigma <- gaussian.mle$cov
	rm(gaussian.mle)

	## Again compute the EMVE scale based on the Gaussian MLE
	if( is.matrix(Sigma) ){
		#if( rcond(Sigma) >= curMinRcond ){
			a <- .emve.scale.missing.Rcpp( Sigma, miss.group.unique, miss.group.counts )$a
	      	Sigma <- a*Sigma
			
			x_mu_diff <- sweep(x.filled, 2, mu, "-")
		    partial <- .partial.mahalanobis.Rcpp(x_mu_diff, Sigma, miss.group.unique, miss.group.counts)
			scale0 <- .fast.scale.mve(partial, scale0_init$cc, scale0_init$ck)
			
			curMaxScaleInd <- which.max( ss0 )
			if( scale0 < max(ss0, na.rm=T) ){
				ss0[ curMaxScaleInd ] <- scale0
				c0[ curMaxScaleInd, ] <- mu
				C0[ , , curMaxScaleInd] <- Sigma
				curMinRcond <- min(curMinRcond, 10^floor( log10( rcond(Sigma) ) ) )
			}				
		#}
	}

	## Final concentration
	for(ii in 1:n.candidates){
	
		x_mu_diff <- sweep(x.filled, 2, c0[ii,], "-")
		partial <- .partial.mahalanobis.Rcpp(x_mu_diff, C0[,,ii], miss.group.unique, miss.group.counts)
		scale0 <- .fast.scale.mve(partial, scale0_init$cc, scale0_init$ck)
      	partial.corrected <- pchisq(partial/scale0,pp)

      	res <- rank(partial.corrected)
      	mm <- (1:n)[res<=floor(n/2)]
      	x.half <- x[mm,]
  
		gaussian.mle <- .cov.EM.rough(x.half, maxits)
		##gaussian.mle <- cov.EM(x.half)
		mu <- gaussian.mle$mean
		Sigma <- gaussian.mle$cov
		rm(gaussian.mle)

  		if( is.matrix(Sigma) ){
			if( rcond(Sigma) >= curMinRcond ){
				a <- .emve.scale.missing.Rcpp( Sigma, miss.group.unique, miss.group.counts )$a
				Sigma <- a*Sigma

				x_mu_diff <- sweep(x.filled, 2, mu, "-")
				partial <- .partial.mahalanobis.Rcpp(x_mu_diff, Sigma, miss.group.unique, miss.group.counts)
				scale0 <- .fast.scale.mve(partial, scale0_init$cc, scale0_init$ck)
				if(scale0 < ss0[ii]){
					ss0[ii] <- scale0
					c0[ii,] <- mu
					C0[,,ii] <- Sigma
				}

				index0 <- order(ss0)
				ss0 <- ss0[index0]
				c0 <- c0[index0,]
				C0 <- C0[,,index0]
			}
		}
	}
	
	## Obtain the best 5 based on MVE scale
	rcondition <- rep(0, n.candidates)

	for(ii in 1:n.candidates){ 
      	C0[,,ii] <- ss0[ii]*C0[,,ii]
      	rcondition[ii] <- rcond(C0[,,ii]) 
	}

	#index0 <- order(rcondition,decreasing = T)

	## newly added
	ss0 <- ss0[1:nf]
	c0 <- c0[1:nf,]
	C0 <- C0[,,1:nf]
	
	## compute pmd
	x.tmp <- x - matrix( c0[1,], n, p, byrow=T)
	pmd <- .partial.mahalanobis.Rcpp(x.tmp, C0[,,1], miss.group.unique, miss.group.counts)
	
	## need to reorder x, pmd for output
	id.reorder <- order(id.order)
	x <- x[id.reorder,]
	pmd <- pmd[ id.reorder]
	pu <- pp[ id.reorder ] 
	
	## compute adjusted pmd
	pmd.adj <- qchisq( pchisq( pmd, df=pu, log.p=T, lower.tail=F), df=p, log.p=T, lower.tail=F) 
	pmd.adj[ which( pu == p)] <- pmd[ which(pu==p) ]

	res <- list(
		S = C0[,,1],
		mu = c0[1,],
		mve.scale=ss0[1],
		cand.mve.scale=ss0,
		cand.mu=c0,
		cand.Sigma=C0,
		pmd = pmd,
		pmd.adj = pmd.adj,
		pu = pu)
	res	
}


