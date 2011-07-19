\name{BaB}
\alias{BaB}

\title{
	Finding the exact p-value.
}
\description{
	\code{BaB} finds an exact p-value by solving a 0-1 knapsack problem.
	The 0-1 knapsack problem is solved by a branch and bound algorithm.
	For more details, see \cite{Higgins, Rivest, Stark}.
}
\usage{
	BaB(Z, t = NULL, asTaint = FALSE, asNumber = FALSE, 
		M = NULL, takeOutZeroMMB=TRUE, give.strategy = FALSE,
		bound.col = "e.max", calc.e_p=calc.pairwise.e_p, 
		 w_p = weight.function("no.weight"))
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
	\item{give.strategy}{If \code{give.strategy = TRUE}, output will include the solution
		to the 0-1 knapsack problem.  
	}
	 \item{bound.col, calc.e_p, w_p}{Arguments used to compute 
	 	\code{t} from audit data, instead of passing \code{t} directly.  
		These arguments are ignored if \code{t} is not NULL.
		See \code{\link{compute.stark.t}} for details. }
}
\details{
	\code{BaB} pre-processes the data to make the branch and bound algorithm
	more efficient, and obtains all information from \code{Z} necessary to perform
	the branch and bound algorithm.
	\code{BaB} then calls \code{\link{runBaB}}, which calls the branch and bound 
	function.
	
	When \code{give.strategy = TRUE}, the output of the solution will be 
	a vector \code{strategy} of size \code{length(nrow(Z$strat))}.
	The solution can be obtained by, for each stratum \code{i}, putting
	\code{e.max} amount of difference in the \code{strategy[i]} batches
	corresponding to the largest values of \code{u}.  
	For more details, see \cite{Higgins, Rivest, Stark}.
}
\seealso{
  See \code{\link{LKPBound}} 
  for finding a p-value through a continuous relaxation.  
  See \code{\link{eqValBound}} and \code{\link{withReplaceBound}} 
  for finding a p-value through other relaxations.  
  See \code{\link{runBaB}} for running the branch and bound algorithm given
  a value vector \code{u}, a cost vector \code{q}, a margin \code{M}, and
  a \code{CIDnum} vector.
  See \code{\link{compute.stark.t}} for computing \code{t} through audit data.
}

\examples{
	data(MN_Senate_2006)
	BaB(MN_Senate_2006.strat, takeOutZeroMMB = FALSE, give.strategy = TRUE)
}	

\author{
	Mike Higgins, Hua Yang
}

\references{
  M. Higgins, R. L. Rivest, P. B. Stark. 
  \emph{Sharper p-Values for Stratified Election Audits} 
}
