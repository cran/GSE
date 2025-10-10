
generate.randcorr <- function(cond, p, tol=1e-5, maxits=100) {
	if( p < 3 ) stop("p must be larger than 2.")
	lambda <- sort(c(1, runif(p-2, min=1, max=cond), cond))
	x <- matrix(rnorm(p*p),p,p)
	Sigma <- x%*%t(x)
	Q <- eigen(Sigma, symmetric=TRUE)$vectors
	ratio <- 0
	iter <- 0
	while (abs(ratio - cond) > tol & iter < maxits) {
		iter <- iter + 1
		Sigma <- Q%*%diag(lambda)%*%t(Q)
		Sigma <- diag(diag(Sigma)^(-1/2))%*%Sigma%*%diag(diag(Sigma)^(-1/2))
		eS <- eigen(Sigma, symmetric=TRUE)
		Q <- eS$vectors
		lambda <- eS$values
		ratio <- lambda[1]/lambda[p]
		lambda[p] <- lambda[1]/cond
	}
	return(Sigma)
}

.generate.clean <- function(n, p, cond, A=NULL){
	if( is.null(A) ) A <- generate.randcorr(cond, p)
	x <- mvrnorm(n, mu=rep(0, p), Sigma=A)
	return(list(x=x, A=A))
}

generate.cellcontam <- function(n, p, cond, contam.size, contam.prop, A=NULL){
	x <- .generate.clean(n, p, cond, A)
	contam.num <- floor(n*p*contam.prop)
	u <- matrix( 0, n, p)
	if( contam.num > 0){
		u[sample(1:(n*p), contam.num)] <- 1
		x$x[ which(u == 1)] <- contam.size + rnorm(contam.num, sd=0.01)
	}
	x$u <- u
	x
}

generate.casecontam <- function(n, p, cond, contam.size, contam.prop, A=NULL){
	x <- .generate.clean(n, p, cond, A)
	Aeig <- eigen(x$A, symmetric=T)$vector
	Aevec <- Aeig[,p]
	Aevec.size <- sqrt(drop(Aevec%*%solve(x$A)%*%Aevec))
	Aevec <- contam.size*Aevec/Aevec.size
	contam.num <- floor(n*contam.prop)
	u <- matrix(0, n, p)
	if( contam.num > 0){
		##u[ sample(1:n, contam.num), ] <- 1
		u[ 1:contam.num, ] <- 1
		x$x[ rowSums(u) == p] <- matrix(Aevec, contam.num, p, byrow=T) + rnorm(contam.num*p, sd=0.01)
	}
	x$u <- u
	x
}

generate.mixedcontam <- function(n, p, cond, contam.size, contam.prop, A=NULL){
  x <- .generate.clean(n, p, cond, A)
  Aeig <- eigen(x$A, symmetric=T)$vector
  Aevec <- Aeig[,p]
  Aevec.size <- drop(sqrt(Aevec%*%solve(x$A)%*%Aevec))
  Aevec <- contam.size[1]*Aevec/Aevec.size
  contam.num1 <- floor(n*contam.prop[1])
  u <- matrix(0, n, p)
  if( contam.num1 > 0){
    ##u[ sample(1:n, contam.num), ] <- 1
    u[ 1:contam.num1, ] <- 1
    x$x[ rowSums(u) == p] <- matrix(Aevec, contam.num1, p, byrow=T) + rnorm(contam.num1*p, sd=0.01)
  }
  n2 <- n*(1-contam.prop[1])
  y <- x$x[(n*contam.prop[1] +1):n,]
  contam.num2 <- floor(n2*p*contam.prop[2])
  u2 <- matrix( 0, n2, p)
  if( contam.num2 > 0){
    u2[sample(1:(n2*p), contam.num2)] <- 1
    y[ which(u2 == 1)] <- contam.size[2] + rnorm(contam.num2, sd=0.01)
  }
  x$u <- rbind(u[1:contam.num1,], u2)
  x$x[(n*contam.prop[1] +1):n,] <- y
  x
}

