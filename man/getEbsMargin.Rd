\name{getEbsMargin}
\alias{getEbsMargin}

\title{
	Updating Error Bounds and Margin Given the Observed Maximum
}
\description{
	Gives updated values \code{u} and \code{M} 
	given a value of the observed maximum \code{t}.
}
\usage{
	getEbsMargin(Z, t, asTaint = FALSE, asNumber = FALSE, M = NULL)
}
\arguments{
	\item{Z}{ A \code{\link{strat.elec.data}} object.}
	\item{t}{ Value of the observed maximum, either as the MRO, as taint, or as the overstatement of the margin in votes.}
	\item{asTaint}{ Set \code{asTaint = TRUE} if \code{t} is the maximum observed taint.}
	\item{asNumber}{ Set \code{asNumber} if \code{t} is the maximum observed overstatement of the margin in votes. }
	\item{M}{ \emph{A priori} margin. If NULL, \code{M} defaults to 1.}
}
\details{
	Creates values \code{u} and margin \code{M} that can 
	be passed into the branch and bound function.  
	The following definitions for \code{u} and \code{M} are described 
	in \cite{Higgins, Rivest, Stark}. 
	The quantity \code{e.max} is obtained through \code{\link{maximumMarginBound}}.
	\itemize{
		\item{Default}{ \code{u = e.max - min(e.max, t). M =  M - sum(min(e.max, t)).}}
		\item{\code{asTaint = TRUE}}{\code{ u = e.max*(1 - t). M = M - sum(e.max*t)}}
		\item{\code{asNumber = TRUE}}{ Same as Default with \code{ t  = t/Z$Margin}.}	}
	
	The output of \code{getEbsMargin} is a list consisting of
	\itemize{
		\item{\code{M}}{The updated margin.}	
		\item{\code{u}}{The updated value vector.}
	}
}

\examples{
	data(MN_Senate_2006)
	getEbsMargin(MN_Senate_2006.strat, t = 0.009, asTaint = TRUE)
}	

\author{
	Mike Higgins, Hua Yang
}

\references{
  M. Higgins, R. L. Rivest, P. B. Stark. 
  \emph{Sharper p-Values for Stratified Election Audits} 
}
