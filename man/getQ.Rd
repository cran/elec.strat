\name{getQ}
\alias{getQ}

\title{
	Obtaining the Cost Vector.
}
\description{
	Obtains the cost value \code{q}, which can then be passed into the branch and bound function.  See \cite{Higgins, Rivest, Stark} for details.
}
\usage{
	getQ(Z)
}
\arguments{
	\item{Z}{ A \code{\link{strat.elec.data}} object.}
}
\examples{
	data(MN_Senate_2006)
	getQ(MN_Senate_2006.strat)
}	

\author{
	Mike Higgins, Hua Yang
}

\references{
  M. Higgins, R. L. Rivest, P. B. Stark. 
  \emph{Sharper p-Values for Stratified Election Audits} 
}
