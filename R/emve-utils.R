###################################################################
## Compute Gaussian MLE of mean and covariance of a dataset by 
## calling cov.EM written by Mike in his robmiss package. This function
## allows users to specify the maximum iteration when computing the MLE
.cov.EM.rough <- function(x, maxits){
    unitw <- rep(1,nrow(x))
    prep.x <- .prelim.norm(x, w1=unitw, w2=unitw)
    estim.em <- NULL
    try(estim.em <- .em.norm(prep.x, showits=FALSE, maxits=maxits),silent=TRUE)
    if (!is.null(estim.em))
      {
        eme <- .getparam.norm(prep.x,estim.em)
      }
    else
      eme <- NULL
    res <- list()
    res$mean <- eme$mu
    res$cov <- eme$s
    return(res)
}

###################################################################
## Compute the correction constants k and 50% chi-sq quantiles for 
## the computation of the scale of MVE for incomplete data (new version)
.scale.mve.init <- function(pp){
	cc <- qchisq(1/2,pp)
	k1 <- 1/pp
	k2 <- 0.5^(pp/2)
	k3 <- 1/(gamma(pp/2))
	k4 <- cc^(1+pp/2)
	k5 <- exp(-cc/2)
	kk <- k1 * k2 * k3 * k4 * k5
	ck <- cc*kk
	return( list(cc=cc, ck=ck) )
}

###################################################################
## Compute the scale of MVE for incomplete data (new version)
## Note: this new version is faster than the old as it only requires
## to compute both the correction constants k and 50% chi-sq quantiles once
.fast.scale.mve <- function(d,cc,ck)
{
	ee <- d/cc
	ff <- order(ee)
	ee <- ee[ff]
	ck <- ck[ff]
	cct <- sum(ck)
	ck <- ck/cct
	ckcumsum <- cumsum(ck)
	n <- length(d)
	tt <- (ckcumsum >= 0.5)
	gg <- 1:n
	rr <- gg[tt]
	mm <- min(rr)
	if( mm==1){
		me <- ee[1]
	} else if(ckcumsum[mm]==0.5){
		me <- ee[mm]
	} else{
		aa <- 1/(0.5-ckcumsum[mm-1])
		bb<-1/(ckcumsum[mm]-0.5)
		dd<-aa+bb
		aa<-aa/dd
		bb<-bb/dd
		me<-aa*ee[mm-1]+bb*ee[mm]
	}
	me
}

###################################################################
## Perform resampling and allow the subsamples to have its covariance matrix 
## satisfied a certain minimum recipricol condition number 
.emve.resamp <- function(x, x.nonmiss, n.resample, n.sub.size, minRcond){
	res <- tryCatch( .Call("emve_resamp", x, x.nonmiss, n.resample, n.sub.size, minRcond, package="GSE"),
		"std::range_error" = function(e){
		conditionMessage( e ) } )
	if( is.character( res ) )
		stop( paste(res, "\nMay consider increase subsample size. The current subsample size is: ", n.sub.size, sep="") )
	return(res)
}

###################################################################
## Compute EMVE scale 
.emve.scale.missing.Rcpp <- function(Sigma, miss_group_unique, miss_group_counts){
	.Call("emve_scale_missing", Sigma, miss_group_unique, miss_group_counts, package="GSE")
}

###################################################################
## Partial mahalanobis distance
## Rcpp version
.partial.mahalanobis.Rcpp <- function(x_mu_diff, Sigma, miss_group_unique, miss_group_counts){
	res <- tryCatch( .Call("fast_partial_mahalanobis", x_mu_diff, Sigma, miss_group_unique, miss_group_counts),
		"std::range_error" = function(e){
		conditionMessage( e ) } )
	if( is.character( res ) ) stop(res)
	return( c(res) )
}






	
	
	