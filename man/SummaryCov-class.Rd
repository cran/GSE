\name{SummaryCov-class}
\docType{class}
\alias{SummaryCov-class}
\alias{show,SummaryCov-method}

\title{Class "SummaryCov" - displaying summary of "CovRobMiss" objects }
\description{ Displays summary information for \code{\link{CovRobMiss-class}} objects }
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("SummaryCov", ...)}.
}
\section{Slots}{
    \describe{
        \item{\code{obj}:}{\code{\link{CovRobMiss-class}} object}
        \item{\code{evals}:}{Eigenvalues and eigenvectors of the covariance or correlation matrix}
    }
}
\section{Methods}{ 
  \describe{
    \item{show}{\code{signature(object = "SummaryCov")}: display the object }
  }
}
\author{ Andy Leung \email{andy.leung@stat.ubc.ca} }
\keyword{classes}

