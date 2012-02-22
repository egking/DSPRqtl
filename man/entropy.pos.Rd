\name{entropy.pos}
\alias{entropy.pos}
\title{Entropy at a position}
\usage{
  entropy.pos(peakChr, peakPos, phenotype.dat, design)
}
\arguments{
  \item{peakChr}{character vector of length one. Must be
  one of the major chromosome arms in the \emph{Drosophila}
  genome ('X','2L','2R','3L',or '3R').}

  \item{peakPos}{numeric vector of length one. A position
  in base pairs in the DSPR position list (every 10kb).}

  \item{phenotype.dat}{\code{data.frame} containing a
  column of ril ids (must be named patRIL) and phenotypes.}

  \item{design}{a character string. One of either 'inbredA'
  or 'inbredB' corresponding to the pA and pB set of inbred
  RILs. Other crossing designs will be supported in the
  future.}
}
\value{
  A numeric vector of length one: the entropy at the given
  position for the set of RILs in the phenotype.dat
  \code{data.frame}.
}
\description{
  \code{entropy.pos} calculates the entropy (proportion of
  missing information) at a given position.
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}
\references{
  Shannon, C.E. 1984. A mathematical theory of
  communication. \emph{Bell System Technical Journal}
  27(3): 379-423.
  \url{http://en.wikipedia.org/wiki/Information_theory#Entropy}
}

