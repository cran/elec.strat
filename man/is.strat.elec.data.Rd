\name{is.strat.elec.data}
\alias{is.strat.elec.data}

\title{
	Verifying a strat.elec.data Object.
}
\description{
	Verifies that an object is a \code{\link{strat.elec.data}} object.
}
\usage{
	is.strat.elec.data(Z)
}
\arguments{
	\item{Z}{ An object.  \code{is.strat.elec.data} is TRUE when \code{Z} is a \code{\link{strat.elec.data}} object.}
}
\seealso{
  See \code{\link{strat.elec.data}} or \code{\link{makeStratObj}} for building a \code{\link{strat.elec.data}} object.  
}
\examples{
	data(MN_Senate_2006)
	is.strat.elec.data(MN_Senate_2006.strat)
}	

\author{
	Mike Higgins
}
