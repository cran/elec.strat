\name{get.first.r.samp}
\alias{get.first.r.samp}
\alias{get.next.r.samp}
\alias{get.prop.samp}

\title{
	Obtain a Vector of Sample Sizes Given Constraint on p-Value
}	
\description{
	\code{get.first.r.samp, get.next.r.samp}, and \code{get.prop.samp} 
	obtain sample sizes so that, if a maximum observed overstatement
	of \code{t} or less is observed, the sample will produce a p-value
	less than \code{alpha}.
	
	\code{get.first.r.samp} uses the \code{\link{first.r}} algorithm to obtain the sample,
	 \code{get.next.r.samp} uses the \code{\link{next.r}} algorithm to obtain the sample, 
	and \code{get.prop.samp} finds a vector of sample sizes that is proportional
	to stratum sizes.
		
	For details about the \code{\link{first.r}} and the
	\code{\link{next.r}} algorithms, and for a description on how
	to produce a sample that will ensure that the p-value 
	is less than \code{alpha}
	when no overstatement greater than \code{t} is uncovered, 
	see \cite{Higgins, Rivest, Stark}.
}
\usage{
	get.first.r.samp(Z, alpha, t, bal=TRUE, numSamp = TRUE, initn = 1, 
		asTaint = FALSE, asNumber = FALSE, M = NULL, 
		takeOutZeroMMB=TRUE)
	get.next.r.samp(Z, alpha, t, bal=TRUE, numSamp = TRUE, initn = 1, 
		asTaint = FALSE, asNumber = FALSE, M = NULL, 
		takeOutZeroMMB=TRUE)
	get.prop.samp(Z, alpha, t, bal=TRUE, numSamp = TRUE, initn = 1, 
		asTaint = FALSE, asNumber = FALSE, M = NULL, 
		takeOutZeroMMB=TRUE)
}
\arguments{
	\item{Z}{ A \code{\link{strat.elec.data}} object.}
	\item{t}{ Value of the observed maximum, either as the MRO, 
		as taint, or as the overstatement of the margin in votes.}
	\item{alpha}{ Threshold for the p-value.  
		If an audit does not uncover an overstatement less than \code{t},
		the sample obtained will ensure that the p-value is less than \code{alpha}.}
	\item{bal}{ If \code{bal = TRUE}, the output will include the 
		expected number of audited ballots for the sample.}
	\item{numSamp}{ If \code{numSamp = TRUE}, the output will include the 
		total number of audited batches.}
	\item{initn}{The first sample size checked by algorithm will have
		 a total of \code{initn} samples. 
		 If this first sample will not produce a p-value less than alpha, 
		 the algorithm will increment the number of samples
		 until such a vector of sample sizes is found.
		 \code{initn} may be adjusted to dramatically decrease the runtime of
		 algorithms.}
	\item{asTaint}{ Set \code{asTaint = TRUE} if \code{t} is the maximum observed taint.}
	\item{asNumber}{ Set \code{asNumber} if \code{t} 
		is the maximum observed overstatement of the margin in votes. }
	\item{M}{ \emph{A priori} margin. If NULL, \code{M} defaults to 1.}
	\item{takeOutZeroMMB}{ 
		Setting \code{takeOutZeroMMB = TRUE} will consider
		 batches with a \code{\link{maximumMarginBound}} of zero as having no 
		 chance of being sampled.}
}
\details{	
	Sample sizes from \code{get.first.r.samp} and \code{get.next.r.samp}
	are obtained by repeatedly calling \code{\link{first.r}} and \code{\link{next.r}},
	respectively,
	while incrementing the total number of samples \code{n}. 
	The algorithm stops when the sample produced will ensure a 
	p-value less than \code{alpha}.
}
\seealso{
  See \code{\link{first.r}, \link{next.r}} 
  and \code{\link{propSizes}} for finding sample sizes 
  given constraints on the p-value and the largest observed overstatement. 
  Also, see \code{\link{first.r}} and \code{\link{next.r}} for a brief 
  description of the \code{\link{first.r}} and \code{\link{next.r}} algorithms.
  See \code{\link{optStrat}} for finding optimal sample sizes so that, 
  if a maximum observed overstatement
	of \code{t} or less is observed, the sample will produce a p-value
	less than \code{alpha}.
 Optimal sample sizes will minimize the number of batches required for audit.

}

\examples{
	data(CA_House_2008)
	get.first.r.samp(CA_House_2008.strat[[3]], alpha = .1, t = .01, asTaint = TRUE)
	get.next.r.samp(CA_House_2008.strat[[3]], alpha = .1, t = .01, asTaint = TRUE)
	get.prop.samp(CA_House_2008.strat[[3]], alpha = .1, t = .01, asTaint = TRUE)
}	

\author{
	Mike Higgins
}

\references{
  M. Higgins, R. L. Rivest, P. B. Stark. 
  \emph{Sharper p-Values for Stratified Election Audits} 
}
