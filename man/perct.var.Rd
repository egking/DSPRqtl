\name{perct.var}
\alias{perct.var}
\title{Effect size}
\usage{
  perct.var(peakChr, peakPos, model, phenotype.dat, design)
}
\arguments{
  \item{peakChr}{character vector of length one. Must be
  one of the major chromosome arms in the \emph{Drosophila}
  genome ('X','2L','2R','3L',or '3R').}

  \item{peakPos}{numeric vector of length one. A position
  in base pairs in the DSPR position list (every 10kb).}

  \item{model}{an object of class formula: a symbolic
  description of the null model to be fitted at each
  position (e.g., \code{phenotype ~ 1}).  The genotype
  effects to be fitted will be added based on
  \code{design}.}

  \item{phenotype.dat}{\code{data.frame} containing a
  column of ril ids (must be named patRIL) and phenotypes.}

  \item{design}{a character string. One of either 'inbredA'
  or 'inbredB' corresponding to the pA and pB set of inbred
  RILs. Other crossing designs will be supported in the
  future.}
}
\value{
  A numeric vector of length one: the percent variance
  explained by the QTL.
}
\description{
  Effect size
}
\details{
  \code{perct.var} calculates the percent of the variance
  explained by the QTL.
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

