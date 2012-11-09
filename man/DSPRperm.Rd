\name{DSPRperm}
\alias{DSPRperm}
\alias{print.pt}
\title{DSPR permutation test}
\usage{
  DSPRperm(model, design, phenotype.dat, batch = 1000,
    niter = 1000, alpha = 0.05)
}
\arguments{
  \item{model}{an object of class formula: a symbolic
  description of the null model to be fitted at each
  position (e.g., \code{phenotype \~ 1}). The genotype
  effects to be fitted will be added based on
  \code{design}.}

  \item{design}{a character string. One of either 'inbredA'
  or 'inbredB' corresponding to the pA and pB set of inbred
  RILs.  Other crossing designs will be supported in the
  future.}

  \item{phenotype.dat}{\code{data.frame} containing a
  column of ril ids (must be named patRIL) and phenotypes.}

  \item{batch}{A numeric vector of length one specifying
  the number of positions to be examined at a time. A
  larger number will use more memory but can be faster.
  Default is 1000.}

  \item{niter}{A numeric vector of length one specifying
  the number of permutations to run. Default is 1000.}

  \item{alpha}{The alpha level for the genome-wide
  signficance threshold. Default is 0.05. Raw maximum LOD
  scores are also returned and \code{\link{quantile}} can
  be used to test other alphas.}
}
\value{
  A list of class \code{pt} containing: \item{maxLODs}{A
  vector containing the maximum LOD score obtained for each
  permutation of the data.} \item{alpha}{the specified
  alpha level} \item{threshold}{the signficance threshold
  at the specified alpha level}
}
\description{
  \code{DSPRperm} performs a permutation test for a DSPR
  dataset.
}
\details{
  The permutation test can take a very long time to run
  (over 24 hrs).
}
\author{
  Elizabeth King (\email{egking@uci.edu})
}

