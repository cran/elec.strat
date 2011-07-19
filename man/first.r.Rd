\name{first.r}
\alias{first.r}
\alias{next.r}
\alias{propSizes}

\title{
	Obtain a Vector of Sample Sizes with Total Number of Samples Fixed
}
\description{
	\code{first.r, next.r}, and \code{propSizes} obtain sample sizes 
	so that the total number of samples is fixed. 
	\code{first.r} uses the \code{first.r} algorithm, 
	\code{next.r} uses the \code{next.r} algorithm, 
	and \code{propSizes} finds a vector of sample sizes that is proportional
	to stratum sizes.
		
	See \cite{Higgins, Rivest, Stark} for details about the \code{first.r} and the
	\code{next.r} algorithms.
}
\usage{
	first.r(Z, n, t = 0, asTaint = FALSE, asNumber = FALSE, M = NULL, initSamp = NULL)
	next.r(Z, n, t = 0, asTaint = FALSE, asNumber = FALSE, M = NULL, initSamp = NULL)
	propSizes(Z, n)
}
\arguments{
	\item{Z}{ A \code{\link{strat.elec.data}} object.}
	\item{n}{ The fixed number of samples. 
		When \code{initSamp} is provided, \code{first.r} and \code{next.r} 
		will run for \code{n} iterations, adding samples iteratively to \code{initSamp};
		\code{first.r} and \code{next.r} will produce a vector of sample sizes 
		with a total of \code{sum(initSamp) + n)} samples.}
	\item{t}{ Value of the observed maximum, either as the MRO, 
		as taint, or as the overstatement of the margin in votes.}
	\item{asTaint}{ Set \code{asTaint = TRUE} if \code{t} is the maximum observed taint.}
	\item{asNumber}{ Set \code{asNumber} if \code{t} is the maximum observed overstatement of the margin in votes. }
	\item{M}{ \emph{A priori} margin. If NULL, \code{M} defaults to 1.}
	\item{initSamp}{ An initial choice of sample sizes.  
		Used in call of \code{\link{get.first.r.samp}} and 
		\code{\link{get.next.r.samp}} to reduce computational time.}
}
\details{
	The arguments \code{t, asTaint, asNumber, M} are used in \code{first.r} and
	\code{next.r} in the call of \code{\link{getEbsMargin}}.  
	The \code{\link{getQ}} function is bypassed to increase efficiency.
	
	\code{propStrat} obtains a vector of sample sizes that has exactly \code{n} samples.
	It obtains such a sample by sorting values of
	\code{k*sum(Z$strat$n)/Z$strat$n}, where \code{k = 0, 1, ...}, in increasing order 
	and allocating a sample to the strata corresponding to the first \code{n} values.  
	Ties are broken by choosing the strata with the largest number of batches.
	See \cite{Higgins, Rivest, Stark} for details.
}
\seealso{
  See \code{\link{get.first.r.samp}, \link{get.next.r.samp}}, 
  and \code{\link{get.prop.samp}} for finding sample sizes 
  given constraints on the p-value and the largest observed overstatement. 
}

\examples{
	data(MN_Senate_2006)
	MN_Senate_2006.strat$strat$audit <- first.r(MN_Senate_2006.strat, n = 150)
	BaB(MN_Senate_2006.strat)
	MN_Senate_2006.strat$strat$audit  <- next.r(MN_Senate_2006.strat, n = 150)
	BaB(MN_Senate_2006.strat)
	MN_Senate_2006.strat$strat$audit  <- propSizes(MN_Senate_2006.strat, n = 150)
	BaB(MN_Senate_2006.strat)
}	

\author{
	Mike Higgins, Hua Yang
}

\references{
  M. Higgins, R. L. Rivest, P. B. Stark. 
  \emph{Sharper p-Values for Stratified Election Audits} 
}
