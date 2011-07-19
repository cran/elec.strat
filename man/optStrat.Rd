\name{optStrat}
\alias{optStrat}

\title{
	Obtain an Optimal Vector of Sample Sizes Given Constraint on p-Value
}	
\description{
	\code{optStrat} will obtain sample sizes so that, 
	if a maximum observed overstatement
	of \code{t} or less is observed, the sample will produce a p-value
	less than \code{alpha}.
	The sample that \code{optStrat} obtains minimizes the total number of
	batches required for audit.
	\code{optStrat} includes options so that, given the number of samples 
	required for audit for optimal sample sizes, 
	the sample that minimizes the expected number of audited ballots is found.
	
	\code{optStrat} can be a very computationally expensive function, and
	should only be used for small contests.
}
\usage{
	optStrat(Z,alpha, t,  bal=TRUE, optBal=FALSE, numSamp = TRUE, 
	asTaint = FALSE, asNumber = FALSE, M = NULL, takeOutZeroMMB=TRUE)
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
	\item{optBal}{ If \code{bal = TRUE}, given the number of batches
		required for audit in an optimal sample, \code{optSamp} will
		find the sample that minimizes the expected number of audited ballots.
		This may dramatically increase the runtime of \code{optStrat}.}
	\item{numSamp}{ If \code{numSamp = TRUE}, the output will include the 
		total number of audited batches.}
	\item{asTaint}{ Set \code{asTaint = TRUE} if \code{t} is the maximum observed taint.}
	\item{asNumber}{ Set \code{asNumber} if \code{t} 
		is the maximum observed overstatement of the margin in votes. }
	\item{M}{ \emph{A priori} margin. If NULL, \code{M} defaults to 1.}
	\item{takeOutZeroMMB}{ 
		Setting \code{takeOutZeroMMB = TRUE} will consider
		 batches with a \code{\link{maximumMarginBound}} of zero as having no 
		 chance of being sampled.}
}

\seealso{
	See \code{\link{get.first.r.samp}, \link{get.next.r.samp}}, 
	and \code{\link{get.prop.samp} }
	for other methods to obtain sample sizes so that, 
	if a maximum observed overstatement
	of \code{t} or less is observed, the sample will produce a p-value
	less than \code{alpha}.
	\code{\link{get.first.r.samp}} uses the \code{\link{first.r}} 
	algorithm to obtain the sample,
	 \code{\link{get.next.r.samp}} uses the \code{\link{next.r}} 
	 algorithm to obtain the sample, 
	and \code{\link{get.prop.samp}} finds a vector of sample sizes that is proportional
	to stratum sizes.
		
}

\examples{
	data(CA_House_2008)
	optStrat(CA_House_2008.strat[[3]], alpha = .1, t = .01, asTaint = TRUE)
	optStrat(CA_House_2008.strat[[3]], alpha = .1, t = .01,
		asTaint = TRUE, optBal = TRUE)
}	

\author{
	Mike Higgins
}

