% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perct.var.R
\name{perct.var}
\alias{perct.var}
\title{Effect size}
\usage{
perct.var(peakChr, peakPos, model, design, phenotype.dat, id.col, sex)
}
\arguments{
\item{peakChr}{character vector of length one. Must be one of the 
major chromosome arms in the \emph{Drosophila} genome 
('X','2L','2R','3L',or '3R').}

\item{peakPos}{numeric vector of length one. A position in base 
pairs in the DSPR position list (every 10kb).}

\item{model}{an object of class formula: a symbolic description of
the null model to be fitted at each position (e.g., 
\code{phenotype ~ 1}).  The genotype effects to be fitted will 
be added based on \code{design}.}

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

\item{sex}{a character string (either 'm' or 'f') specifying the 
sex of the measured individuals. This argument must be supplied 
for a cross design for correct specification of the genotypes on
the X chromosome.}
}
\value{
A numeric vector of length one: the percent variance 
  explained by the QTL.
}
\description{
\code{perct.var} calculates the percent of the variance explained
by the QTL.
}
\author{
Elizabeth King (\email{kingeg@missouri.edu})
}
