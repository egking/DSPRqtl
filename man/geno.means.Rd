\name{geno.means}
\alias{geno.means}
\title{Founder genotype means and standard errors}
\usage{
  geno.means(peakChr, peakPos, model, phenotype.dat,
    design)
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
  A \code{data.frame} of the estimated mean and standard
  error for each founder genotype. If a covariate is
  included in the model statement, the estimate will be the
  founder genotype mean after correcting for the covariate.
}
\description{
  \code{geno.means} estimates the founder genotype means
  and standard errors at a given position.
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

