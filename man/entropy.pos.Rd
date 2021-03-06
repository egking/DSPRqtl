% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy.pos.R
\name{entropy.pos}
\alias{entropy.pos}
\title{Entropy at a position}
\usage{
entropy.pos(peakChr, peakPos, design, phenotype.dat, id.col, sex)
}
\arguments{
\item{peakChr}{character vector of length one. Must be one of the 
major chromosome arms in the \emph{Drosophila} genome 
('X','2L','2R','3L',or '3R').}

\item{peakPos}{numeric vector of length one. A position in base 
pairs in the DSPR position list (every 10kb).}

\item{design}{a character string. For inbred RIL designs:
'inbredA', 'inbredB'. For cross designs: AAcross, BBcross, or
'ABcross'. A and B refer to the pA and pB set of DSPR RILs.}

\item{phenotype.dat}{\code{data.frame} containing phenotype data. 
For inbred designs, there must be a column of numeric RIL ids 
(must be named patRIL). For the ABcross design, there must be
both a patRIL and matRIL column specifying the pA and pB RIL
ids. Cross designs also require a sex column for correct 
specification of the genotypes on the X chromosome.}

\item{id.col}{a character string identifying the name of the 
column containing unique ids for the samples. e.g. for an inbred
design, the patRIL column can be used as the id.}

\item{sex}{a character string (either 'm' or 'f') specifying the 
sex of the measured individuals. This argument must be supplied 
for the AB cross design for correct specification of the
genotypes on the X chromosome.}
}
\value{
A numeric vector: the entropy at the given position for
  the set of RILs in the phenotype.dat \code{data.frame}. In the
  case of the ABcross, the A and B entropy are calculated
  separately and both are returned.
}
\description{
\code{entropy.pos} calculates the entropy (proportion of missing 
information) at a given position.
}
\references{
Shannon, C.E. 1984. A mathematical theory of
  communication. \emph{Bell System Technical Journal} 27(3):
  379-423. 
  \url{http://en.wikipedia.org/wiki/Information_theory#Entropy}
}
\author{
Elizabeth King (\email{kingeg@missouri.edu})
}
