\name{DSPRpeaks}
\alias{DSPRpeaks}
\alias{print.peaks}
\title{Find and Summarize QTL}
\usage{
  DSPRpeaks(qtldat, threshold = 6.8, LODdrop = 2)
}
\arguments{
  \item{qtldat}{An object of class gscan. Output from
  \code{\link{DSPRscan}}.}

  \item{threshold}{numeric vector of length one consisting
  of the signficance threshold. Default is 6.8. Use
  \code{\link{DSPRperm}} to get a threshold specific to a
  given dataset.}

  \item{LODdrop}{numeric vector of length one consisting of
  the LOD drop to be used. Default value is 2 which
  approximates a 95\% confidence interval.}
}
\value{
  A list of class \code{peaks} containing a list for each
  significant peak, each containing: \item{threshold}{the
  specified signficance threshold} \item{peak}{A single row
  \code{data.frame} with the chromosome, physical position
  (bp), genetic position (cM) and LOD score for the peak}
  \item{LODdrop}{the specified LOD drop for the support
  interval} \item{CI}{the upper and lower bounds of the
  confidence interval} \item{founderNs}{the number of RILs
  with each founder genotype at the peak}
  \item{geno.means}{the estimated means and standard errors
  for each founder genotype} \item{perct.var}{the percent
  variation explained by the QTL} \item{entropy}{the
  proporiton of missing information (entropy)}
}
\description{
  \code{DSPRpeaks} takes output from
  \code{\link{DSPRscan}}. Locates and summarizes QTL peaks.
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

