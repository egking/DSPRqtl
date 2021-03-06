% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DSPRscan.R
\name{DSPRscan}
\alias{DSPRscan}
\title{DSPR Genome Scan}
\usage{
DSPRscan(model, design, phenotype.dat, id.col, batch = 1000, sex)
}
\arguments{
\item{model}{an object of class formula: a symbolic description of
the null model to be fitted at each position (e.g., 
\code{phenotype ~ 1}). The genotype effects to be fitted will be
added based on \code{design}.}

\item{design}{a character string. One of either 'inbredA', 
'inbredB', or 'ABcross' corresponding to the pA and pB set of 
inbred RILs or the pA-pB cross design. For round robin designs 
or other cross designs, use the more flexible DSPRgenos and 
standard model fitting functions in R.}

\item{phenotype.dat}{\code{data.frame} containing phenotype data. 
For inbred designs, there must be a column of numeric RIL ids 
(must be named patRIL). For the ABcross design, there must be
both a patRIL and matRIL column specifying the pA and pB RIL
ids.}

\item{id.col}{a character string identifying the name of the 
column containing unique ids for the samples. e.g. for an inbred
design, the patRIL column can be used as the id.}

\item{batch}{A numeric vector of length one specifying the number 
of positions to be examined at a time. A larger number will use 
more memory but can be faster. Default is 1000.}

\item{sex}{a character string (either 'm' or 'f') specifying the 
sex of the measured individuals. This argument must be supplied 
for a cross design for correct specification of the genotypes on
the X chromosome.}
}
\value{
A list of class \code{gscan} containing:
\item{LODscores}{A \code{data.frame} consisting of the chromosome,
physical position (bp), genetic position (cM) and LOD score for 
each position.}
\item{model}{the model specification}
\item{design}{the design specification}
\item{phenotype}{the phenotype \code{data.frame} specified}
}
\description{
Function to perform a genome scan for data generated with the DSPR
RILs. This function is capable of fitting standard models. Users 
requiring more flexibility should use the \code{\link{DSPRgenos}} 
function along with standard model fitting functions in R.
}
\author{
Elizabeth King (\email{kingeg@missouri.edu})
}
