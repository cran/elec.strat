\name{runBaB}
\alias{runBaB}

\title{
	Calling the Branch and Bound Algorithm
}
\description{
	\code{runBaB} calls the branch and bound algorithm.
	The branch and bound algorithm is coded in \code{C}.
}
\usage{
	runBaB(u,q,M,CIDnum)
}
\arguments{
	\item{u}{ A vector of values.  Can be obtained through \code{\link{getEbsMargin}}.}
	\item{q}{ A vector of costs.  Can be obtained through \code{\link{getQ}}.}
	\item{M}{The constraint on the values.  
		Can be obtained through \code{\link{getEbsMargin}}.}
	\item{CIDnum}{ A vector that gives the CIDnum identification for each batch.  
		Can be found at \code{Z$V$CIDnum}. }
}
\examples{
	data(MN_Senate_2006)
	M.u <- getEbsMargin(MN_Senate_2006.strat, t=2, asNumber = TRUE)
	u <- M.u$u
	M <- M.u$M
	q <- getQ(MN_Senate_2006.strat)
	CIDnum <- MN_Senate_2006.strat$V$CIDnum
	runBaB(u,q,M,CIDnum)
}	

\author{
	Mike Higgins
}

