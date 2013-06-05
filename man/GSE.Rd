\name{GSE}
\alias{GSE}
\title{ Generalized S-Estimator in the presence of missing data}
\description{
  Computes the Generalized S-Estimate (GSE) -- a robust estimate of location and scatter 
  for data with contamination and missingness. 
}
\usage{
GSE(x, tol=1e-5, maxiter=500, init="emve", tol.scale=1e-4, miter.scale=30, 
	print.step=1, mu0, S0, ...)
}
\arguments{
  \item{x}{a matrix or data frame. May contain missing values, but cannot contain columns with completely missing entries. }
  \item{tol}{tolerance for the convergence criterion. Default is 1e-5.}
  \item{maxiter}{maximum number of iterations for the GSE algorithm. Default is 500.}
  \item{init}{
  type of initial estimator. Currently this can either be "emve" (emve), "qc", or "huber" 
  (see \code{\link{emve}} and \code{\link{HuberPairwise}}). Default is "emve". 
  If \code{mu0} and \code{S0} are provided, this argument is ignored.}
  \item{tol.scale}{tolerance for the computation of the GS-scale. Default is 1e-4.}
  \item{miter.scale}{maximum number of iterations for the computation the GS-scale. Default is 30.}
  \item{print.step}{this argument determines the level of printing which is done during the estimation process. 
	The default value is 1 means that the number of iteration at convergence is printed out;
	a value of 0 means no printing occurs, except for error messages; and a value of 2 means that the 
	relative change of scale and estimated GSE scale at each iteration is printed out.}
  \item{mu0}{optional vector of initial location estimate}
  \item{S0}{optional matrix of initial scatter estimate}
  \item{...}{optional arguments for computing the initial estimates (see \code{\link{emve}}, \code{\link{HuberPairwise}}).}
}
\details{
   This function computes GSE as described in Danilov et al. (2012). The estimator requires a robust positive definite 
   initial estimator. This initial estimator is required to ``re-scale" the partial square mahalanobis distance for 
   the different missing pattern, in which a single scale parameter is not enough. This function currently allows two 
   main different initial estimators: EMVE (the default; see \code{\link{emve}} and Huberized Pairwise 
   (see \code{\link{HuberPairwise}}). GSE using the latter estimator when the tuning constant
   \eqn{c0} is 0 is referred to as QGSE in Danilov et al. (2012). Numerical results have shown that GSE with EMVE as 
   initial has better performance (in both efficiency and robustness), but computing time can be longer. 
}
\value{
  An S4 object of class \code{\link{GSE-class}} which is a subclass of the virtual class \code{\link{CovRobMissSc-class}}. The
  output S4 object contains the following slots:
  
	\tabular{ll}{
	\code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
	\code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
	\code{sc} \tab Generalized S-scale (GS-scale). Can be accessed via \code{\link{getScale}}. \cr
	\code{pmd} \tab Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. \cr
	\code{pmd.adj} \tab Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. \cr
	\code{pu} \tab Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. \cr
	\code{iter} \tab Number of iterations till convergence. Not meant to be accessed. \cr
	\code{eps} \tab relative change of the GS-scale at convergence. Not meant to be accessed. \cr
	\code{call} \tab Object of class \code{"language"}. Not meant to be accessed. \cr
	\code{x} \tab Input data matrix. Not meant to be accessed. \cr
	\code{p} \tab Column dimension of input data matrix. Not meant to be accessed. \cr
	\code{estimator}\tab Character string of the name of the estimator used. Not meant to be accessed. \cr
	}
}
\references{
 	Danilov, M., Yohai, V.J., Zamar, R.H. (2012). Robust Esimation of Multivariate
	Location and Scatter in the Presence of Missing Data. \emph{Journal of the American Statistical Association}
	\bold{107}, 1178--1186.
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Ruben H. Zamar, Mike Danilov, Victor J. Yohai}
\seealso{
  \code{\link{emve}}, \code{\link{HuberPairwise}}, \code{\link{GSE-class}}
}
\examples{
\dontrun{
set.seed(12345)
n <- 100; p <- 10
A <- matrix(0.9, p, p); diag(A) <- 1
x <- mvrnorm(n, rep(0,p), A)

## Introduce 10\% contamination points 
pcont <- 0.1; dcont <- 10
A.svd <- svd(A)
U <- A.svd$u \%*\% diag( sqrt(A.svd$d) ) \%*\% t(A.svd$v)

## Contaminated points in the direction corresponding to smallest eigenvalue
vv <- A.svd$u[,p]
uu <- 1/sqrt( t(vv) \%*\% solve(A, vv) )
vv <- vv*uu
ncont <- rbinom(n,1,pcont)
if( sum(ncont) > 0 )
	x[which(ncont== 1),] <- dcont * matrix(vv,sum(ncont),p,byrow=T)

## Introduce 10\% missingness
pmiss <- 0.1
nmiss <- matrix(rbinom(n*p,1,pmiss), n,p)
x[ which( nmiss == 1 ) ] <- NA

## Using EMVE as initial
res.emve <- GSE(x)
summary(res.emve) ## summary of the output 
plot(res.emve) ## plot of the output
slrt( getScatter(res.emve), A) ## LRT distances to the true covariance

## Using QC as initial
res.qc <- GSE(x, init="qc")
summary(res.qc)
plot(res.qc)
slrt( getScatter(res.qc), A) ## in general performs worse than if EMVE used as initials
}
}

