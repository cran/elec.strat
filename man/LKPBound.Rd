\name{LKPBound}
\alias{LKPBound}

\title{
	p-value Through a Continuous Relaxation
}
\description{
	Finds a p-value through the LKP Bound: a continuous relaxation bound of the 
	original 0-1 knapsack problem.  
	Offers an option to include a lower-bound in output, thus 
	computing an upper and lower bound on the exact p-value.
	See \cite{Higgins, Rivest, Stark} for more details.
}
\usage{
	LKPBound(Z, t = NULL, asTaint = FALSE, asNumber = FALSE, 
		M = NULL, takeOutZeroMMB = TRUE, LKP.lower.bound = FALSE, 
		bound.col = "e.max",  
		calc.e_p=calc.pairwise.e_p, w_p = weight.function("no.weight"))	
}
\arguments{
	\item{Z}{ A \code{\link{strat.elec.data}} object.}
	\item{t}{ Value of the observed maximum, either as the MRO, 
		as taint, or as the overstatement of the margin in votes.}
	\item{asTaint}{ Set \code{asTaint = TRUE} if \code{t} is the maximum observed taint.}
	\item{asNumber}{ Set \code{asNumber} if \code{t} is the maximum observed overstatement of the margin in votes. }
	\item{M}{ \emph{A priori} margin. If NULL, \code{M} defaults to 1.}
	\item{takeOutZeroMMB}{ 
		Setting \code{takeOutZeroMMB = TRUE} will consider
		 batches with a \code{\link{maximumMarginBound}} of zero as having no 
		 chance of being sampled.}
	 \item{LKP.lower.bound}{ 
	 	Set \code{LKP.lower.bound = TRUE} to compute 
	 	a lower-bound of the exact p-value in addition to the upper-bound computed
		from the continuous relaxation.  
	 	Lower-bound computed according to \cite{Higgins, Rivest, Stark}.}
	 \item{bound.col, calc.e_p, w_p}{Arguments used to compute 
	 	\code{t} from audit data, instead of passing \code{t} directly.  
		These arguments are ignored if \code{t} is not NULL.
		See \code{\link{compute.stark.t}} for details. }
}
\seealso{
  See \code{\link{eqValBound}} and \code{\link{withReplaceBound}} 
  for finding a p-value through other relaxations.  
  See \code{\link{BaB}} for finding an exact p-value 
  through solving a 0-1 knapsack problem. 
  See \code{\link{compute.stark.t}} for computing \code{t} through audit data.
}

\examples{
	data(MN_Senate_2006)
	LKPBound(MN_Senate_2006.strat, takeOutZeroMMB = FALSE)
	LKPBound(MN_Senate_2006.strat, t = 2, asNumber = TRUE, 
		takeOutZeroMMB = FALSE, LKP.lower.bound = TRUE)
	
	data(CA_House_2008)	
	CA_House_2008.strat[[1]]$strat$audit <- 1
	LKPBound(CA_House_2008.strat[[1]], t = 0, LKP.lower.bound = TRUE)
}	

\author{
	Mike Higgins, Hua Yang
}

\references{
  M. Higgins, R. L. Rivest, P. B. Stark. 
  \emph{Sharper p-Values for Stratified Election Audits} 
}
