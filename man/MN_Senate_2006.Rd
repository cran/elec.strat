\name{MN_Senate_2006}
\alias{MN_Senate_2006}
\alias{MN_Senate_2006.strat}

\docType{data}
\title{
  2006 Minnesota U.S. Senate Race
}
\description{
  Contains a \code{\link{strat.elec.data}} object for the 
  2006 Minnesota U.S. Senate Race named \code{MN_Senate_2006.strat}.
  
}
\usage{data(MN_Senate_2006)}
\details{
  The winner of the election was Amy Klobuchar. 
  Mark Kennedy was the
  runner-up. 
  There were a total of 2,217,818 voters, and the margin of
  victor was 443,196 votes. 
  The largest precinct wise difference between
  the hand count and machine count was a 2-vote swing from Amy Klobuchar
  to Mark Kennedy. 
}
\references{
  M. Halvorson and L. Wolff. \emph{Report and analysis of the 2006
    post-election audit of Minnesotas voting systems}.
  \url{http://ceimn.org/files/CEIMNAuditReport2006.pdf}
}
\examples{
data(MN_Senate_2006)
BaB(MN_Senate_2006.strat)
}
\keyword{datasets}
