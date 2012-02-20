\name{DSPRscan}
\alias{DSPRscan}
\title{DSPR Genome Scan}
\usage{
  DSPRscan(model, design, phenotype.dat, batch = 1000)
}
\arguments{
  \item{model}{an object of class formula: a symbolic
  description of the null model to be fitted at each
  position (e.g., \code{phenotype ~ 1}). The genotype
  effects to be fitted will be added based on
  \code{design}.}

  \item{design}{a character string. One of either 'inbredA'
  or 'inbredB' corresponding to the pA and pB set of inbred
  RILs. Other crossing designs will be supported in the
  future.}

  \item{phenotype.dat}{\code{data.frame} containing a
  column of ril ids (must be named patRIL) and phenotypes.}

  \item{batch}{A numeric vector of length one specifying
  the number of positions to be examined at a time. A
  larger number will use more memory but can be faster.
  Default is 1000.}
}
\value{
  A list of class \code{gscan} containing:
  \item{LODscores}{A \code{data.frame} consisting of the
  chromosome, physical position (bp), genetic position (cM)
  and LOD score for each significant peak.}
  \item{model}{the model specification} \item{design}{the
  design specification} \item{phenotype}{the phenotype
  \code{data.frame} specified}
}
\description{
  Function to perform a genome scan for data generated with
  the DSPR RILs.
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

