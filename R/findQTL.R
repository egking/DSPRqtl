##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param qtldat 
##' @param threshold 
##' @return 
##' @author Elizabeth King
findQTL<-function(qtldat,threshold)
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
