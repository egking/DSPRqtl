\name{findCI}
\alias{findCI}
\title{Calculate LOD support interval}
\usage{
  findCI(peakChr, peakPos, qtldat, peak, LODdrop)
}
\arguments{
  \item{peakChr}{character vector of length one. Must be
  one of the major chromosome arms in the \emph{Drosophila}
  genome ('X','2L','2R','3L',or '3R').}

  \item{peakPos}{numeric vector of length one. A position
  in base pairs in the DSPR position list (every 10kb).}

  \item{qtldat}{\code{data.frame} of chromosome, position,
  and LOD scores (column names chr,Ppos,Gpos,LOD). List
  element LODscores from \code{\link{DSPRscan}}.}

  \item{peak}{numeric vector of length one consisting of
  the peak LOD score.}

  \item{LODdrop}{numeric vector of length one consisting of
  the LOD drop to be used. Default value is 2 which
  approximates a 95\% confidence interval.}
}
\value{
  A \code{data.frame} with two rows containing the
  chromosome, physical position (bp), genetic position (cM)
  and LOD scores corresponding to the lower and upper
  bound.
}
\description{
  \code{findCI} calculates the LOD support interval for a
  given LOD peak.
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

