\name{makeStratObj}
\alias{makeStratObj}

\title{
	Making a strat.elec.data Object from an elec.data Object
}
\description{
	Makes a \code{\link{strat.elec.data}} object from an \code{\link{elec.data}} object.
}
\usage{
	makeStratObj(Z, strat.col = NULL, CID = NULL, auditTable = NULL)
}
\arguments{
	\item{Z}{ An \code{\link{elec.data}} object.}
	\item{strat.col}{ Name of column in \code{Z$V} that identifies strata.  If no value of strat.col is passed, \code{makeStratObj} will assume that stratum ID is contained in \code{Z$V$CID}.}
	\item{CID}{ A vector of length \code{nrow(Z$V)} that identifies strata.}
	\item{auditTable}{ A \code{data.frame} of dimension \code{length(unique(CID)} x 2 used to identify the number of samples taken from each stratum. Including \code{auditTable} is not necessary if \code{Z} contains audit information.}
}
\details{
	\code{makeStratObj} requires as input a valid \code{\link{elec.data}} object 
	\code{Z} such that one of the following is true:
	\itemize{
		\item \code{Z} has a vector \code{Z$V$CID} that identifies strata.
		\item A \code{strat.col} name identifying the stratification column in \code{Z$V} is passed to \code{makeStratObj}.
		\item A \code{CID} vector of length \code{nrow(Z$V)} specifying the stratification is passed to \code{makeStratObj}.
	}
	If \code{Z$audit} is \code{NULL}, information giving the number 
	of sampled batches in each stratum can be included through \code{auditTable}.
	The argument \code{auditTable} should be a 
	\code{data.frame} of dimensions \code{unique(CID)} x 2.  	
	The first column is a list of unique stratum IDs.
	The second column is the number of batches sampled 
	within the corresponding stratum.
	
	\code{makeStratObj} creates a \code{data.frame} \code{Z$strat}.
	The columns of \code{Z$strat} are
	\itemize{
		\item{ \code{CID}}{ The ID of the stratum.}
		\item{ \code{CIDnum}}{ A number between 1 and \code{length(Z$strat$CID)} assigned to that stratum.}
		\item{ \code{n}}{ The number of batches contained in that stratum.}
		\item{\code{ audit}}{ The number of batches sampled from that stratum.}
	}
	If no audit data is provided, \code{Z$strat$audit} defaults to a zero vector.

	If \code{Z$V$CID} is NULL, \code{makeStratObj} will copy the stratum 
	labels into \code{Z$V$CID}.  
	\code{makeStratObj} will also create
	\itemize{
		\item{\code{Z$CID.col}}{\code{"CID"}}
		\item{\code{Z$CIDnum}}{ The \code{CIDnum} of the stratum.}
		\item{\code{Z$V$e.max}}{ \code{maximumMarginBound(Z)}}
	}
	If \code{Z$audit[Z$PID.col]} is not NULL, \code{makeStratObj} will create
	\code{Z$audit$e.max}, the \code{maximumMarginBound(Z)} 
	for batches in \code{Z$audit}.
	
	After sending an \code{\link{elec.data}} object through \code{makeStratObj}, 
	the object will be both
	an \code{\link{elec.data}} object and a \code{\link{strat.elec.data}} object.
}
\seealso{
  See \code{\link{strat.elec.data}} to create a \code{\link{strat.elec.data}} object from a votes \code{data.frame} and an audit \code{data.frame}.
}
\examples{
	data("CA_House_2008")
	dstrat <- CA_House_2008.strat[[1]]
	auditTable <- cbind(unique(dstrat$V$CID),1)
	dstrat <- makeStratObj(dstrat,auditTable = auditTable)
}	

\author{
	Mike Higgins, Hua Yang
}
