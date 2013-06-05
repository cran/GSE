###################################################################
## This file contains functions extracted from robmiss packages
## written by Mike Danilov
slrt <- function(S, trueS){
    ## Standardized LRT-statistic for testing if cov S is the same as trueS
    if (missing(trueS)) trueS <- diag(rep(1,ncol(S)))
    if (is.character(trueS)) trueS <- .std.cov(p=ncol(S),trueS)
    if (!is.matrix(trueS)) stop("True covariance (trueS) misspecified")

    return(.lrt(solve(trueS)%*%S))
}
  
.lrt <- function(S)
{
  ## LRT-statistic for testing if covariance is identity  
  if (sum(is.na(S))>0) return(NA)
  return(sum(diag(S))-log(det(S))-ncol(S))
}
  
.std.cov <- function(p, ct, seed=2727)
{
  if (!exists(".Random.seed")) runif(1) # to initialize .Random.seed
  seed.before <- .Random.seed
  if (!is.null(seed)) set.seed(seed)
  res <- switch(ct,
                h=cor(mvrnorm(p+1, mu=rep(0,p), Sigma=diag(rep(1,p)))),
                m=cor(mvrnorm(round(p+5), mu=rep(0,p), Sigma=diag(rep(1,p)))),
                l=cor(mvrnorm(p*10, mu=rep(0,p), Sigma=diag(rep(1,p)))),
                i=diag(rep(1,p)))
  .Random.seed <- seed.before
  return(res)
}
#***********************************************************************
# Changes NA's to single precision missing value code
.na.to.snglcode <- function(x,mvcode){
      x[is.na(x)] <- as.double(mvcode)
  x}
#***********************************************************************
# Changes missing value code to NA
.code.to.na <- function(x,mvcode){
      x[x==mvcode] <- NA
  x}
#***********************************************************************
#  Perform preliminary manipulations on matrix of continuous data.  
#  Rows are sorted by missing data pattern.
.prelim.norm <- function(x,w1,w2,do.std=TRUE){
  if (missing(w1)) w1 <- rep(1,nrow(x))
  if (missing(w2)) w2 <- rep(1,nrow(x))
# get dimensions of x
  if(is.vector(x)) x <- matrix(x,length(x),1)
  n <- nrow(x); p <- ncol(x); storage.mode(x) <- "double"
  w1n <- as.double(sum(w1)/n)
  storage.mode(w1) <- "double"
  w2n <- as.double(sum(w2)/n)
  storage.mode(w2) <- "double"
  
# find missingness patterns
  r <- 1*is.na(x)
  nmis <- as.integer(apply(r,2,sum))
  names(nmis) <- dimnames(x)[[2]]
  
# index the missing data patterns
  ## MD: this is likely suboptimal but I'm just trying to make as few
  ## changes as possible.
  mdp <- .enum.miss(r) # as.numeric((r%*%(2^((1:ncol(x))-1)))+1)
# do row sort
  ro <- order(mdp)
  x <- matrix(x[ro,],n,p)
  w1 <- w1[ro]
  w2 <- w2[ro]
  mdp <- mdp[ro]
  r <- matrix(r[ro,],n,p)
  ro <- order(ro)

# compress missing data patterns
  mdpst <- as.integer(seq(along=mdp)[!duplicated(mdp)])
  mdp <- unique(mdp); npatt <- length(mdpst)

## create r-matrix for display purposes
  r <- 1-r; r <- matrix(r[mdpst,],npatt,p)
  if(npatt==1) tmp <- format(n)
  if(npatt>1)  tmp <- format(c(mdpst[2:npatt],n+1)-mdpst)
  dimnames(r) <- list(tmp,dimnames(x)[[2]])
  storage.mode(r) <- "integer"

  ## center and scale the columns of x
  if (do.std)
    {
      if(sum(is.na(x))<length(x)){
        mvcode <- as.double(max(x[!is.na(x)])+1000)
        x <- .na.to.snglcode(x,mvcode)
        tmp <- .Fortran("wctrsc",x,w1,n,p,numeric(p),numeric(p),mvcode)
        x <- tmp[[1]]; xbar <- tmp[[5]]; sdv <- tmp[[6]]
        x <- .code.to.na(x,mvcode)}  
      if(sum(is.na(x))==length(x)){
        xbar <- rep(0,p); sdv <- rep(1,p)}
    } else {
      xbar <- rep(0,p); sdv <- rep(1,p) }
  
# form matrix of packed storage indices
  d <- as.integer((2+3*p+p^2)/2)
  psi <- .Fortran("mkpsi",p,matrix(as.integer(0),p+1,p+1))[[2]]
# other bookkeeping quantities
  if(npatt>1) nmdp <- as.integer(c(mdpst[-1],n+1)-mdpst)
  if(npatt==1) nmdp <- n
  sj <- .Fortran("sjn",p,npatt,r,integer(p))[[4]]
  nmon <- .Fortran("nmons",p,npatt,nmdp,sj,integer(p))[[5]]
  last <- .Fortran("lasts",p,npatt,sj,integer(npatt))[[4]]
  tmp <- .Fortran("layers",p,sj,integer(p),integer(1))
  layer <- tmp[[3]]; nlayer <- tmp[[4]]
# return list
  list(x=x,w1=w1,n=n,w1n=w1n,w2=w2,w2n=w2n, p=p,r=r,nmis=nmis,ro=ro,mdpst=mdpst,
    nmdp=nmdp,npatt=npatt,xbar=xbar,sdv=sdv,
    d=d,psi=psi,sj=sj,nmon=nmon,last=last,layer=layer,nlayer=nlayer)}
#***********************************************************************
# Retrieves means and covariances from theta. If corr=FALSE , returns
# a list containing a vector of means and a covariance matrix. If
# corr=TRUE , returns a list containing a vector of means, a vector of
# standard deviations, and a correlation matrix.
.getparam.norm <- function(s,theta,corr=FALSE ){
  mu <- theta[s$psi[1,2:(s$p+1)]]*s$sdv + s$xbar
  names(mu) <- dimnames(s$x)[[2]]
  sigma <- theta[s$psi[2:(s$p+1),2:(s$p+1)]]
  sigma <- matrix(sigma,s$p,s$p)
  tmp <- matrix(s$sdv,s$p,s$p)
  sigma <- sigma*tmp*t(tmp)
  dimnames(sigma) <- list(names(mu),names(mu))
  if(corr){
    sdv <- sqrt(diag(sigma)); names(sdv) <- names(mu)
    tmp <- matrix(sdv,s$p,s$p)
    r <- sigma/(tmp*t(tmp)); dimnames(r) <- list(names(mu),names(mu))
    result <- list(mu=mu,sdv=sdv,r=r)}
  else result <- list(mu=mu,sigma=sigma)
  result}
#***********************************************************************
# Makes a theta vector out of a list of specified parameters.
.makeparam.norm <- function(s,thetalist){
  result <- numeric(s$d); result[1] <- -1
  xbar <- s$xbar;sdv <- s$sdv
  mu <- (thetalist[[1]]-xbar)/sdv
  result[2:(s$p+1)] <- mu
  if(length(thetalist)==3){
    tmp <- matrix(thetalist[[2]],s$p,s$p)
    sigma <- thetalist[[3]]*tmp*t(tmp)}
  else sigma <- thetalist[[2]]
  tmp <- matrix(sdv,s$p,s$p)
  sigma <- sigma/(tmp*t(tmp))
  tmp <- as.vector(s$psi[2:(s$p+1),2:(s$p+1)])
  result[tmp] <- as.vector(sigma)
  result}
#***********************************************************************
# Finds posterior mode of theta under the multivariate
# normal model. If no prior is specified, finds the mle.
.em.norm <- function(s,start,showits=TRUE ,maxits=1000,criterion=.0001,
     prior){
  s$x <- .na.to.snglcode(s$x,999)
  if(missing(start)){
    start <- .Fortran("stvaln",s$d,numeric(s$d),s$p,s$psi)[[2]]}
  tmp <- as.integer(numeric(s$p))
  tobs <- .Fortran("wtobsn",s$d,numeric(s$d),s$p,s$psi,s$n,s$w1n,s$w2n,s$x,s$w1,s$npatt,
    s$r,s$mdpst,s$nmdp,tmp)[[2]]
# iterate to mle
  it <- 0; converged <- FALSE
  if(showits) cat(paste("Iterations of EM:","\n"))
  while((!converged)&(it<maxits)){
  old <- start

  start <- .Fortran("wemn",s$d,old,start,tobs,s$p,s$psi,s$n,s$w1n,s$w2n,
    s$x,s$w1,s$w2,s$npatt,s$r,s$mdpst,s$nmdp,tmp,tmp,numeric(s$p))[[3]]
  
# print iteration number
  it <- it+1; if(showits) cat(paste(format(it),"...",sep=""))
#  print(det(getparam.norm(s,start)$sigma))
  converged <- max(abs(old-start))<=criterion}
  if(showits)cat("\n")
  start}


##***********************************************************************

## added by Mike Danilov, 2008-02-20
.enum.miss <- function(mm)
  {
    ord <- do.call("order",as.data.frame(mm))
    enum <- c(0,cumsum(apply(abs(mm[ord[-1],,drop=FALSE]-mm[ord[-nrow(mm)],,drop=FALSE]),1, sum)>0))
    enum[ord] <- enum
    return(enum)
  }

CovEM <- function(x)
  {
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

	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")

	xcall <- match.call()
  
    unitw <- rep(1,nrow(x))
    prep.x <- .prelim.norm(x, w1=unitw, w2=unitw)
    estim.em <- NULL
    try(estim.em <- .em.norm(prep.x, showits=FALSE),silent=TRUE)
    if (!is.null(estim.em))
      {
        eme <- .getparam.norm(prep.x,estim.em)#/.6745
      }
    else
      eme <- NULL

	pmd <- partial.mahalanobis(x.orig, eme$mu, eme$s)
	res <- new("CovRobMiss",
		call = xcall,
		S = eme$s,
		mu = eme$mu,
		estimator = "Maximum likelihood estimator", 
		x = x.orig,
		pmd = pmd@pmd,
		pmd.adj = pmd@pmd.adj,
		p = pmd@p,
		pu = pmd@pu)
	res
	
  }
