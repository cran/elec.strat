\name{CA_House_2008}
\alias{CA_House_2008}
\alias{CA_House_2008.strat}

\docType{data}
\title{
  Set of 2008 California U.S. House Races
}
\description{
  A \code{\link{list}} of 20 \code{\link{strat.elec.data}} objects, each
  containing voting data for a contest in the 2008 California U.S. House Race.
  The data for contest \code{i} is contained in the \code{\link{strat.elec.data}}
  object \code{CA_House_2008.strat[[i]]}.
}
\usage{
	data(CA_House_2008)
}
\details{
	Each of the 20 contests had exactly two candidates that received 
	a large portion of the vote.  
	Each contest was contained within 2 to 5 counties.
	\code{\link{optStrat}} can find sample sizes for most of these contests in
	a reasonable amount of time.
		
	Data for the contests were obtained through the 
	California Statewide Database (SWDB).
	The data can be found at 
	\url{http://swdb.berkeley.edu/pub/data/G08/state/state_g08_sov_data_by_g08_svprec.dbf}.
}
\examples{
data(CA_House_2008)
optStrat(CA_House_2008.strat[[3]], alpha = .1, t = 0)
}
\keyword{datasets}
