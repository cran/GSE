\name{emve-class}
\docType{class}
\alias{emve-class}
\alias{getCandidates,emve-method}
\title{Extended Minimum Volume Ellipsoid (EMVE) in the presence of missing data.}
\description{ Class of Extended Minimum Volume Ellipsoid. It has the superclass of \code{CovRobMissSc}. }
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("emve", ...)},
    but the best way of creating \code{emve} objects is a call to the function
    \code{emve} which serves as a constructor.
}
\section{Slots}{
	\describe{
	\item{\code{mu}}{Estimated location. Can be accessed via \code{\link{getLocation}}. }
	\item{\code{S}}{Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. } 
    \item{\code{sc}}{Estimated EMVE scale. Can be accessed via \code{\link{getScale}}. }
    \item{\code{cand.sc}}{EMVE scales for the top 10 candidates. Can be accessed via \code{\link{getCandidates}}. }
    \item{\code{cand.mu}}{a matrix of dimension 10 by p of the estimated EMVE location for the top 10 candidates, where p is the dimension of the data. Can be accessed via \code{\link{getCandidates}}. }
    \item{\code{cand.S}}{an array of dimension p by p by 10 of the estimated EMVE scatter for the top 10 candidates. Can be accessed via \code{\link{getCandidates}}. }		
	\item{\code{pmd}}{Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. }
	\item{\code{pmd.adj}}{Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. }
	\item{\code{pu}}{Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. }
	\item{\code{call}}{Object of class \code{"language"}. Not meant to be accessed. }
	\item{\code{x}}{Input data matrix. Not meant to be accessed. }
	\item{\code{p}}{Column dimension of input data matrix. Not meant to be accessed. }
	\item{\code{estimator}}{Character string of the name of the estimator used. Not meant to be accessed. }
	}	
}
\section{Extends}{ Class \code{"\linkS4class{CovRobMissSc}"}, directly. }
\section{Methods}{ 
The following methods are defined with the superclass "CovRobMiss":
    \describe{
    \item{show}{\code{signature(object = "CovRobMiss")}: display the object }
    \item{summary}{\code{signature(object = "CovRobMiss")}: calculate summary information }
	\item{plot}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: plot the object. See \code{\link{plot}}}
	\item{getDist}{\code{signature(object = "CovRobMiss")}: return the squared partial Mahalanobis distances }
	\item{getDistAdj}{\code{signature(object = "CovRobMiss")}: return the adjusted squared partial Mahalanobis distances }
	\item{getDim}{\code{signature(object = "CovRobMiss")}: return the dimension of observed entries for each case }
	\item{getLocation}{\code{signature(object = "CovRobMiss")}: return the estimated location vector}
	\item{getScatter}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: return the estimated scatter matrix}
	\item{getMissing}{\code{signature(object = "CovRobMiss")}: return the case number(s) with completely missing data, if any}
	\item{getOutliers}{\code{signature(object = "CovRobMiss", cutoff = "numeric")}: return the case number(s) adjusted squared 
		distances above \code{(1 - cutoff)}th quantile of chi-square p-degrees of freedom.}
	}
	In addition to above, the following methods are defined with the class "CovRobMissSc":
	  \describe{
 	\item{getScale}{\code{signature(object = "CovRobMissSc")}: return the MVE scale of the best candidate}
	}
	In addition to above, the following methods are defined with the class "emve":
	  \describe{
 	\item{getCandidates}{\code{signature(object = "emve")}: return a list of the top 10 candidates of EMVE estimates:
		\code{cand.sc} (a vector of EMVE scales), \code{cand.mu} (a matrix of EMVE location), and \code{cand.S}
		(an array of EMVE scatter). See \code{\link{emve}} for details.}
}
}
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }
\seealso{
   \code{\link{emve}}, \code{\link{CovRobMissSc-class}}, \code{\link{CovRobMiss-class}}
}
\keyword{classes}


