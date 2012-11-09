##' \code{findQTL} finds peak LOD scores that exceed the given
##' threshold (identifies signficant QTL).
##' 
##' @title Find QTL
##'   
##' @param qtldat \code{data.frame} of chromosome, position, and LOD
##'   scores (column names chr,Ppos,Gpos,LOD). List element LODscores 
##'   from \code{\link{DSPRscan}}.
##'   
##' @param threshold numeric vector of length one consisting of the 
##'   signficance threshold. Default is 6.8. Use
##'   \code{\link{DSPRperm}} to get a threshold specific to a given
##'   dataset.
##'   
##'   
##' @return A \code{data.frame} consisting of the chromosome, physical
##'   position (bp), genetic position (cM) and LOD score for each
##'   significant peak.
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##'
findQTL<-function(qtldat,threshold=6.8)
{
  output<-data.frame('chr'=numeric(length=0),'Ppos'=numeric(length=0),'Gpos'=numeric(length=0),'LOD'=numeric(length=0))
  for(i in 1:nrow(qtldat))
  {
    if(i==1 | i==nrow(qtldat)){}else{
      if(qtldat$LOD[i]>qtldat$LOD[i-1] & qtldat$LOD[i]>qtldat$LOD[i+1] & qtldat$LOD[i]>threshold)
      {
        peak<-data.frame('chr'=qtldat[i,'chr'],'Ppos'=qtldat[i,'Ppos'],'Gpos'=qtldat[i,'Gpos'],'LOD'=qtldat[i,'LOD'])
        output<-rbind(output,peak)
      }
    }
  }
  return(output)
}
