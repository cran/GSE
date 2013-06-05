\name{plot-methods}
\docType{methods}
\alias{plot}
\alias{plot-method}
\alias{plot,CovRobMiss-method}
\alias{plot,CovRobMiss,missing-method}
\title{Plot methods for objects of class 'CovRobMiss'}
\description{
    Plot methods for objects of class 'CovRobMiss'. The following plots are available:
	
	- chi-square qqplot for adjusted square partial Mahalanobis distances
	
	- index plot for adjusted square partial Mahalanobis distances

	- distance-distance plot comparing the adjusted distances based on classical MLE and robust estimators
	
	Cases with completely missing data will be dropped out. Outliers are identified using some pre-specific 
	cutoff value, for instance 99\% quantile of chi-square with p degrees of freedom, where p is the column
	dimension of the data. Identified outliers can also be retrieved using \code{\link{getOutliers}} with
	an optional argument of \code{cutoff}, ranged from 0 to 1.
}
\usage{
	\S4method{plot}{CovRobMiss}(x, which = c("all","distance","qqchi2", "dd"), 
		which=c("all", "distance", "qqchisq", "dd"),
		ask = (which=="all" && dev.interactive(TRUE)),
		cutoff, ...)
} 
\arguments{
 \item{x}{an object of class \code{"CovRobMiss"}}
  \item{which}{ Which plot to show? Default is \code{which}="all". }
  \item{ask}{ logical; if 'TRUE', the user is \emph{ask}ed before each plot, see 'par(ask=.)'. 
	Default is \code{ask = which=="all" && dev.interactive()}.  }
  \item{cutoff}{ The cutoff value for the distances.  }
  \item{...}{ Additional arguments to be passed over to control the coordinates. See \code{\link[ggplot2]{coord_trans}}.  }
}
\examples{
\dontrun{
data(boston)
res <- GSE(boston)

## plot all graphs
plot(res)

## plot individuals plots
plot(res, which="qqchisq")
plot(res, which="distance")
plot(res, which="dd")

## control the coordinates, e.g. log10 transform the y-axis
plot(res, which="qqchisq", ytrans="log10", xtrans="log10")
plot(res, which="distance", ytrans="log10")
plot(res, which="dd", ytrans="log10", xtrans="log10")
}}
\keyword{methods}

