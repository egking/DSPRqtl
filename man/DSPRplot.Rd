\name{DSPRplot}
\alias{DSPRplot}
\title{Genome Scan Plot}
\usage{
DSPRplot(qtldata, threshold, legNames = NULL)
}
\arguments{
  \item{qtldata}{a \code{list} of output from DSPRscan with
  each list element corrsponding to one DSPRscan result.
  e.g. to plot multple genome scans on a single plot:
  \code{qtldata <- list(scan1,scan2)}}

  \item{threshold}{numeric vector of length one consisting
  of the signficance threshold. Default is 6.8 for inbred
  designs and 10.1 for the ABcross. Use
  \code{\link{DSPRperm}} to get a threshold specific to a
  given dataset. dataset.}

  \item{legNames}{a character vector with names for each
  DSPRscan result to be plotted. Defaults to the phenotype
  column names used in DSPRscan.}
}
\description{
Function to plot genome scan results for the DSPR RILs
}
\author{
Elizabeth King (\email{egking@uci.edu})
}

