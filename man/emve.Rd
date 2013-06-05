\name{emve}
\alias{emve}
\title{Extended Minimum Volume Ellipsoid (EMVE) in the presence of missing data}
\description{
  Computes the Extended S-Estimate (ESE) version of the minimum volume ellipsoid (EMVE), which is
  used as an initial estimator in Generlized S-Estimator (GSE) for missing data by default. 
}
\usage{emve(x, n.resample=100, maxits=5, n.sub.size) }
\arguments{
  \item{x}{a matrix or data frame. May contain missing values, but cannot contain columns with completely missing entries.}
  \item{n.resample}{integer indicating the number of subsamples. Default is 100. Must be at least 20.}
  \item{maxits}{integer indicating the maximum number of iterations of Gaussian MLE calculation for each subsample. Default is 5.}
  \item{n.sub.size}{optional integer indicating the size of a subsample. Must be at least the column dimension of the data.}
}
\details{
  This function computes EMVE as described in Danilov et al. (2012). Subsampling algorithm is used for computing EMVE.
  By default, we take \eqn{N = 100} subsamples of size \eqn{n_0 = p/(1-\alpha)}, where \eqn{p} is the dimension of the data and
  \eqn{\alpha} is the fraction of missing data. The subsample size \eqn{n_0} must be chosen to be larger than \eqn{p} to avoid singularity. 
  Only the best 10 candidates of estimated location and scatter from \eqn{N} subsamples will be output.
  
  In the algorithm, there exists a concentration step in which Gaussian MLE is computed for \eqn{50\%} of the data points using the 
  classical EM-algorithm multiplied by a scalar factor. This step is repeated for each subsample. As the computation can be heavy as the 
  number of subsample increases, we set by default the maximum number of iteration of classical EM-algorithm (i.e. \code{maxits}) as 5. 
  Users are encouraged to refer to Danilov et al. (2012) for details about the algorithm
  and Rubin and Little (2002) for the classical EM-algorithm for missing data. 
}
\value{
  An S4 object of class \code{\link{emve-class}} which is a subclass of the virtual class \code{\link{CovRobMissSc-class}}.
  The output S4 object contains the following slots:
  
	\tabular{ll}{
	\code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
	\code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
    \code{sc} \tab Estimated EMVE scale. Can be accessed via \code{\link{getScale}}. \cr
    \code{cand.sc} \tab EMVE scales for the top 10 candidates. Can be accessed via \code{\link{getCandidates}}. \cr
    \code{cand.mu} \tab a matrix of dimension 10 by p of the estimated EMVE location for the top 10 candidates, where p is the dimension of the data. Can be accessed via \code{\link{getCandidates}}. \cr
    \code{cand.S} \tab an array of dimension p by p by 10 of the estimated EMVE scatter for the top 10 candidates. Can be accessed via \code{\link{getCandidates}}. \cr		
	\code{pmd} \tab Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. \cr
	\code{pmd.adj} \tab Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. \cr
	\code{pu} \tab Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. \cr
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
	
	Rubin, D.B. and Little, R.J.A. (2002). \emph{Statistical analysis with missing data} (2nd ed.). New York: Wiley.
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Ruben H. Zamar, Mike Danilov, Victor J. Yohai}
\examples{
\dontrun{
set.seed(12345)
n <- 100; p <- 10
A <- matrix(0.9, p, p); diag(A) <- 1
x <- mvrnorm(n, rep(0,p), A)

## Introduce 10% contamination points 
pcont <- 0.1; dcont <- 10
A.svd <- svd(A)
U <- A.svd$u %*% diag( sqrt(A.svd$d) ) %*% t(A.svd$v)

## Contaminated points in the direction corresponding to smallest eigenvalue
vv <- A.svd$u[,p]
uu <- 1/sqrt( t(vv) %*% solve(A, vv) )
vv <- vv*uu
ncont <- rbinom(n,1,pcont)
if( sum(ncont) > 0 )
	x[which(ncont== 1),] <- dcont * matrix(vv,sum(ncont),p,byrow=T)

## Introduce 10% missingness
pmiss <- 0.1
nmiss <- matrix(rbinom(n*p,1,pmiss), n,p)
x[ which( nmiss == 1 ) ] <- NA

## Use EMVE
res <- emve(x)

## returns the list of candidate estimates
res.cand <- getCandidates(res)
str(res.cand)

## the element 'cand.Sigma' as the best 10 scatter
## in arrays with the 3rd dimension corresponds to their
## rankings, e.g. best candidate
res.cand$cand.S[,,1]
res.cand$S ## the best candidate also available as 'mu','S'

## EMVE results could pass on to GSE without recomputing
## the initial estimate as follows:
res.GSE <- GSE(x, mu0=res.cand$mu, S0=res.cand$S)
}
}
\seealso{
  \code{\link{GSE}}, \code{\link{emve-class}}
}
