\name{TSGS}
\alias{TSGS}
\alias{.gy.filt}
\title{Two-Step Generalized S-Estimator for cell- and case-wise outliers}
\description{
  Computes the Two-Step Generalized S-Estimate (2SGS) -- a robust estimate of location and scatter 
  for data with cell-wise and case-wise contamination. 
}
\usage{
TSGS(x, alpha=0.99, it=TRUE, ...)

.gy.filt(x, alpha, it=TRUE)
}
\arguments{
  \item{x}{a matrix or data frame.}
  \item{alpha}{quantile of the reference distribution used in tail comparison in the first step. A standard normal is used as the reference distribution. Default value is 0.99.}
  \item{it}{logical, whether the filtering is repeated until no additional cell-wise outliers are identified in the first step. Default value is TRUE.}
  \item{...}{optional arguments to be used in the computation of GSE in the second step (see \code{\link{GSE}}).}
}
\details{
   This function computes 2SGS as described in Agostinelli et al. (2014). The procedure has two steps:

   In Step I, the method filters (i.e., flags and removes) cell-wise outliers using Gervini-Yohai approach. 
   Outliers are flagged and replaced by missing values (\code{NA}) when the empirical tail distribution is heavier than a reference distribution. 
   A standard normal is currently used as a reference distribution as suggested in Agostinelli et al. (2014). 
   The filtering step can be called on its own by using the function \code{GSE:::.gy.filt(x, alpha, it=FALSE)}.
   
   In Step II, the method applies GSE (see \code{\link{GSE}}), which has been specifically designed to deal with 
   incomplete multivariate data with case-wise outliers, to the filted data coming from Step I.   

   The application to the Chemical data set analyzed in Agostinelli et al. (2014) can be found in \code{\link{geochem}}.   
   
   The tools that were used to generate contaminated data in the simulation study in Agostinelli et al. (2014) can be found in \code{\link{generate.cellcontam}} and \code{\link{generate.casecontam}}.
}
\value{
	The following gives the major slots in the output S4 object: 
	\tabular{ll}{
	\code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
	\code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
	\code{xf} \tab Filtered data matrix from the first step of 2SGS. Can be accessed via \code{\link{getFiltDat}}. \cr
	}
}
\references{
 	Agostinelli, C., Leung, A. , Yohai, V.J., and Zamar, R.H. (2014) Robust estimation of multivariate location and scatter in the presence of cellwise and casewise contamination. arXiv:1406.6031[math.ST]
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Claudio Agostinelli, Ruben H. Zamar, Victor J. Yohai}
\seealso{
  \code{\link{GSE}}, \code{\link{generate.cellcontam}}, \code{\link{generate.casecontam}}
}
\examples{
\dontrun{
set.seed(12345)

# Generate 5\% cell-wise contaminated normal data 
# using a random correlation matrix with condition number 100
x <- generate.cellcontam(n=100, p=10, cond=100, contam.size=5, contam.prop=0.05)

## Using MLE
slrt( cov(x$x), x$A)

## Using Fast-S
slrt( rrcov:::CovSest(x$x)@cov, x$A)

## Using 2SGS
slrt( GSE:::TSGS(x$x)@S, x$A)


# Generate 5\% case-wise contaminated normal data 
# using a random correlation matrix with condition number 100
x <- generate.casecontam(n=100, p=10, cond=100, contam.size=15, contam.prop=0.05)

## Using MLE
slrt( cov(x$x), x$A)

## Using Fast-S
slrt( rrcov:::CovSest(x$x)@cov, x$A)

## Using 2SGS
slrt( GSE:::TSGS(x$x)@S, x$A)

}
}


