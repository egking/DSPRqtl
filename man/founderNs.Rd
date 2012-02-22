\name{founderNs}
\alias{founderNs}
\title{Numbers of founder genotypes}
\usage{
  founderNs(peakChr, peakPos, design, phenotype.dat)
}
\arguments{
  \item{peakChr}{character vector of length one. Must be
  one of the major chromosome arms in the \emph{Drosophila}
  genome ('X','2L','2R','3L',or '3R').}

  \item{peakPos}{numeric vector of length one. A position
  in base pairs in the DSPR position list (every 10kb).}

  \item{design}{a character string. One of either 'inbredA'
  or 'inbredB' corresponding to the pA and pB set of inbred
  RILs. Other crossing designs will be supported in the
  future.}

  \item{phenotype.dat}{\code{data.frame} containing a
  column of ril ids (must be named patRIL) and phenotypes.}
}
\value{
  A named numeric vector consisting of the numbers of RILs
  in the phenotype.dat \code{data.frame} with each founder
  genotype at the given position. A RIl is assigned a
  founder genotype if the HMM probability is over 0.95. The
  number of RILs with a heterozygous genotype and an
  uncertain genotype are also returned.
}
\description{
  \code{founderNs} calculates the number of founder
  genotypes represented in the DSPR RILs at a given
  position.
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

