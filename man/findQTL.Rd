\name{findQTL}
\alias{findQTL}
\title{Find QTL}
\usage{
  findQTL(qtldat, threshold = 6.8)
}
\arguments{
  \item{qtldat}{\code{data.frame} of chromosome, position,
  and LOD scores (column names chr,Ppos,Gpos,LOD). List
  element LODscores from \code{\link{DSPRscan}}.}

  \item{threshold}{numeric vector of length one consisting
  of the signficance threshold. Default is 6.8. Use
  \code{\link{DSPRperm}} to get a threshold specific to a
  given dataset.}
}
\value{
  A \code{data.frame} consisting of the chromosome,
  physical position (bp), genetic position (cM) and LOD
  score for each significant peak.
}
\description{
  \code{findQTL} finds peak LOD scores that exceed the
  given threshold (identifies signficant QTL).
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

