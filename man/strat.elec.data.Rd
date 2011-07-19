\name{strat.elec.data}
\alias{strat.elec.data}

\title{
	Making a strat.elec.data Object from a Votes data.frame and an Audit data.frame
}
\description{
	Makes a \code{strat.elec.data} and an \code{\link{elec.data}} object from 
	a votes \code{data.frame} and an audit \code{data.frame}.
}
\usage{
	strat.elec.data(V, C.names=names(V)[2:length(V)], f = 1, 
	audit=NULL, pool=TRUE, tot.votes.col="tot.votes", PID.col="PID", 
	strat.col = NULL, CID = NULL, auditTable = NULL)
}
\arguments{
	\item{V}{ A \code{data.frame} of votes.}
	\item{C.names}{ Names of candidates.}
	\item{f}{ The number of winners. } 
	\item{audit}{ An audit \code{data.frame}. }
	\item{pool}{ Combine small candidates into single pseudo-candidates to increase power.}
	\item{tot.votes.col}{ Name of column that has the total votes for the batches. }
	\item{PID.col}{ Name of column that identifies unique batches.}
	\item{strat.col}{ Name of column in votes that identifies strata.}  
	\item{CID}{ A vector of length \code{nrow(votes)} that identifies strata.}
	\item{auditTable}{ A \code{data.frame} of dimension \code{length(unique(CID)} x 2 used to identify the number of samples taken from each stratum. The \code{auditTable} is not necessary if an audit \code{data.frame} is included.}
}
\details{
	\code{strat.elec.data} creates a \code{strat.elec.data} object:
	an \code{\link{elec.data}} object with additional entries for easy use with the\code{elec.strat} package.  
	
	\code{strat.elec.data} allows for two ways to specify the stratification:
	\itemize{
		\item Specify \code{strat.col}: the name of the column in \code{V} 
		that contains strata information.
		\item Provide a \code{CID} vector of length \code{nrow(V)} 
		specifying the stratification.
	}
	If neither method is used to specify stratification, only an \code{\link{elec.data}}
	object is created. 

	If \code{audit} is not NULL and \code{strat.col} is provided,
	\code{strat.elec.data} will find the \code{strat.col} column in
	\code{audit} to create an \code{auditTable}. 
	If no column in \code{audit} is labeled as \code{strat.col}, 
	\code{strat.elec.data} will throw an error.

	If \code{audit} is not NULL, \code{strat.col} is NULL, 
	and \code{CID} is provided, 
	\code{strat.elec.data} will try to find 
	the column in \code{audit} labeled \code{PID.col}
	to generate the \code{auditTable}.
	If the \code{PID.col} column is not in \code{audit},
	then \code{auditTable} will need to be given;
	otherwise \code{Z$strat$audit} defaults to a zero vector.
	
	The argument \code{auditTable} should be a 
	\code{data.frame} of dimensions \code{unique(CID)} x 2.  	
	The first column is a list of unique stratum IDs.
	The second column is the number of batches sampled 
	within the corresponding stratum.
	
	\code{strat.elec.data} first calls \code{\link{elec.data}} to create 
	an \code{elec} object \code{Z}.
	An \code{auditTable} is either created through \code{audit}
	or given by the \code{auditTable} argument,
	and \code{\link{makeStratObj}} is called to create
	a \code{strat.elec.data} object. 
	
	For a detailed description of the structure of a \code{strat.elec.data} object,
	see \code{\link{makeStratObj}}.
	For a more detailed description of the arguments 
	\code{V, C.names, f, pool, audit, tot.votes.col}, see \code{\link{elec.data}}.
}
\seealso{
  See \code{\link{elec.data}} to create an \code{\link{elec.data}} object.   
  See \code{\link{makeStratObj}} to create a \code{strat.elec.data} object from an 
   \code{\link{elec.data}} object.   
   Both \code{\link{elec.data}} and \code{\link{makeStratObj}} are 
   called by \code{strat.elec.data}.
}
\examples{
	data(MN_Senate_2006)
	votes <- MN_Senate_2006.strat$V
	audit <- MN_Senate_2006.strat$audit
	CID <- MN_Senate_2006.strat$V$CID
	names <- c("Klo","Ken")
	strat.elec.data(V = votes, C.names = names, audit = audit, CID = CID)
}	

\author{
	Mike Higgins
}
