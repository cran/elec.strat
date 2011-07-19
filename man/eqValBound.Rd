\name{eqValBound}
\alias{eqValBound}
\alias{withReplaceBound}

\title{
	p-value Through Relaxation on Number of Batches Without Difference
}
\description{
	\code{eqValBound} and \code{withReplaceBound} find a p-value by
	changing the original constraint (that the total difference is greater than the margin),
	instead placing a restriction on the
	number of batches with error no larger than \code{t}.
	
	\code{eqValBound} finds an exact solution with this restriction, 
	whereas \code{withReplaceBound} finds a more conservative bound.
	See \cite{Stark} for more details about \code{withReplaceBound}.
}
\usage{
	eqValBound(Z, t = NULL, asTaint = FALSE, asNumber = FALSE, 
		M = NULL, takeOutZeroMMB=TRUE, bound.col = "e.max", 
		 calc.e_p=calc.pairwise.e_p, w_p = weight.function("no.weight"))
	withReplaceBound(Z, t = NULL, asTaint = FALSE, asNumber = FALSE, 
		M = NULL, takeOutZeroMMB = TRUE, bound.col = "e.max",  	
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
	 \item{bound.col, calc.e_p, w_p}{Arguments used to compute 
	 	\code{t} from audit data, instead of passing \code{t} directly.  
		These arguments are ignored if \code{t} is not NULL.
		See \code{\link{compute.stark.t}} for details. }
}
\seealso{
  See \code{\link{LKPBound}} 
  for finding a p-value through a continuous relaxation.  
  See \code{\link{BaB}} for finding an exact p-value 
  through solving a 0-1 knapsack problem. 
  See \code{\link{compute.stark.t}} for computing \code{t} through audit data.
}

\examples{
	data(MN_Senate_2006)
	eqValBound(MN_Senate_2006.strat, takeOutZeroMMB = FALSE)
	withReplaceBound(MN_Senate_2006.strat, t = 2, asNumber = TRUE, 
		takeOutZeroMMB = FALSE)
}	

\author{
	Mike Higgins, Hua Yang
}

\references{
  P.B. Stark. \emph{Conservative Statistical Post-Election Audits}. Annals of Applied Statistics, 2:2. 550-581. 
}
