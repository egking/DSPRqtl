##' DSPRqtl
##' 
##' Analysis of DSPR phenotypes
##' 
##' \tabular{ll}{
##'   Package: \tab DSPRqtl\cr
##'   Type: \tab Package\cr
##'   Version: \tab 2.0-4\cr
##'   Date: \tab 2015-06-21\cr
##'   License: \tab GPL-2\cr
##'   LazyLoad: \tab yes\cr
##'   LazyData: \tab yes\cr
##' }
##' 
##' @name DSPRqtl-package
##' 
##' @aliases DSPRqtl-package DSPRqtl
##' 
##' @docType package
##' 
##' @author Elizabeth King
##' 
##' @keywords package
##' @importFrom graphics abline axis box legend mtext par plot points rect
##' @importFrom stats as.formula dnorm lm logLik optim quantile sd summary.aov
##' @importFrom utils data
##'
NULL

##' Positionlist
##' 
##' List of regularly spaced positions every 10kb used for data 
##' analysis of DSPR data. chr is the chromosome arm, Ppos is the 
##' position in base pairs, Gpos is the position in centiMorgans, and 
##' Gaxis is used for plotting the entire genome on a single axis.
##' 
##' @name positionlist_wgenetic
##' 
##' @aliases positionlist_wgenetic poslist
##' 
##' @docType data
##' 
##' @format A data frame with 4 variables.
##' \describe{
##' \item{\code{chr}}{character vector of chromosome}
##' \item{\code{Ppos}}{numeric vector}
##' \item{\code{Gpos}}{numeric vector}
##' \item{\code{Gaxis}}{numeric vector}
##' }
##' 
##' @keywords datasets
##' 
NULL

##' ADH data
##' 
##' The ADH activity phenotype data.
##' 
##' @name ADHdata
##' 
##' @aliases ADHdata
##' 
##' @docType data
##' 
##' @format A data frame with 3 variables.
##' \describe{
  ##' \item{\code{patRIL}}{numeric vector}
      ##' \item{\code{matRIL}}{numeric vector}
      ##' \item{\code{adh}}{numeric vector}
      ##' }
##' 
##' @keywords datasets
##' 
NULL

##' ADH scan
##' 
##' The ADH genome scan results.
##' 
##' @name ADHscan
##' 
##' @aliases ADHscan
##' 
##' @docType data
##' 
##' @format A list with 4 variables.
##' \describe{
##' \item{\code{LODscores}}{data.frame containing positions and LODscores}
##' \item{\code{model}}{formula}
##' \item{\code{design}}{character vector}
##' \item{\code{phenotype}}{data.frame containing RIL ids and phenotypes}
##' }
##' 
##' @keywords datasets
##' 
NULL
