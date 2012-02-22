##' \code{findCI} calculates the LOD support interval for a given LOD peak.
##'
##'
##' 
##' @title Calculate LOD support interval
##' 
##' @param peakChr character vector of length one. 
##' Must be one of the major chromosome arms in the 
##' \emph{Drosophila} genome ('X','2L','2R','3L',or '3R').
##' 
##' @param peakPos numeric vector of length one. A position in 
##' base pairs in the DSPR position list (every 10kb). 
##'
##' @param qtldat \code{data.frame} of chromosome, position, 
##' and LOD scores (column names chr,Ppos,Gpos,LOD). List element LODscores 
##' from \code{\link{DSPRscan}}.
##' 
##' @param peak numeric vector of length one consisting of the peak LOD score.
##'  
##' @param LODdrop numeric vector of length one consisting of the LOD drop to be used. 
##' Default value is 2 which approximates a 95\% confidence interval.
##' 
##' @return A \code{data.frame} with two rows containing the chromosome, physical position (bp), 
##' genetic position (cM) and LOD scores corresponding to the lower and upper bound. 
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##' 
##' 
findCI<-function(peakChr,peakPos,qtldat,peak,LODdrop)
{
  startIndex<-which(qtldat$chr==peakChr & qtldat$Ppos==peakPos)
  l.ci<-peak-LODdrop  #drop threshold
  start1<-startIndex #this one is for the other side of the peak
  #one side
  while(qtldat$LOD[start1]>l.ci & (start1-1) >= 1) 
  {  
    start1<-start1-1 
  }
  
  #other side
  while(qtldat$LOD[startIndex]>l.ci & (startIndex+1) <= nrow(qtldat))
  {  
    startIndex<-startIndex+1 
  }
  CIs<-qtldat[c(start1,startIndex),]
  rownames(CIs)<-c('Lower Bound','Upper Bound')
  return(CIs)
}
