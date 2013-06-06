##' \code{findCI} calculates the LOD support interval for a given LOD
##' peak.
##' 
##' @title Calculate LOD support interval
##'   
##' @param peakChr character vector of length one. Must be one of the
##'   major chromosome arms in the \emph{Drosophila} genome
##'   ('X','2L','2R','3L',or '3R').
##'   
##' @param peakPos numeric vector of length one. A position in base
##'   pairs in the DSPR position list (every 10kb).
##'   
##' @param qtldat \code{data.frame} of chromosome, position, and LOD
##'   scores (column names chr,Ppos,Gpos,LOD). List element LODscores 
##'   from \code{\link{DSPRscan}}.
##'   
##' @param method a character string specifying the method to use for 
##'   the confidence interval. Options are: 'LODdrop' calculates a LOD 
##'   drop interval for the drop amount specified, 'BCI' calculates the 
##'   Bayesian credible interval for the fraction specified, and 'both' 
##'   calculates and returns both. 
##'   
##' @param LODdrop numeric vector of length one consisting of the LOD
##'   drop to be used when using method 'LODdrop'. Default value is 2 
##'   which approximates a 95\% confidence interval for the inbred 
##'   designs. Users of the ABcross design should consider using a 
##'   larger LOD drop.
##'   
##' @param BCIprob numeric vector of length one consisting of the nominal
##' Bayes fraction. Default value is 0.95. 
##'   
##' @return A \code{data.frame} with two rows containing the
##'   chromosome, physical position (bp), genetic position (cM) and
##'   LOD scores corresponding to the lower and upper bound. A list of
##'   two \code{data.frames} are returned when the method used is 'both'.
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##' 
findCI<-function(peakChr,peakPos,qtldat,method,LODdrop,BCIprob)
{
  if(missing(LODdrop))
  {
    if(method=='LODdrop'|method=='both')
    {
      LODdrop<-2
    }else{
      LODdrop<-NA
    } 
  }  
  
  if(missing(BCIprob))
  {
    if(method=='BCI'|method=='both')
    {
      BCIprob<-0.95
    }else{
      BCIprob<-NA
    } 
  }  
  
  if(method=='LODdrop')
  {
    startIndex<-which(qtldat$chr==peakChr & qtldat$Ppos==peakPos)
    l.ci<-qtldat[startIndex,'LOD']-LODdrop  #drop threshold
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
  }else{
    if(method=='BCI')
    {
      CIdat<-qtldat
      CIdat$chr[CIdat$chr=='2R'|CIdat$chr=='2L']<-'2'
      CIdat$chr[CIdat$chr=='3R'|CIdat$chr=='3L']<-'3'
      
      Bchr<-peakChr
      Bchr[Bchr=='2R'|Bchr=='2L']<-'2'
      Bchr[Bchr=='3R'|Bchr=='3L']<-'3'
      startIndex<-which(qtldat$chr==peakChr & qtldat$Ppos==peakPos)
      BCI<-bayesint(qtldat,Bchr,startIndex,BCIprob)
      rownames(BCI)<-c('Lower Bound','Upper Bound')      
      return(BCI)
      
    }else{
      
      if(method=='both')
      {
        startIndex<-which(qtldat$chr==peakChr & qtldat$Ppos==peakPos)
        l.ci<-qtldat[startIndex,'LOD']-LODdrop  #drop threshold
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
        
        CIdat<-qtldat
        CIdat$chr[CIdat$chr=='2R'|CIdat$chr=='2L']<-'2'
        CIdat$chr[CIdat$chr=='3R'|CIdat$chr=='3L']<-'3'
        
        Bchr<-peakChr
        Bchr[Bchr=='2R'|Bchr=='2L']<-'2'
        Bchr[Bchr=='3R'|Bchr=='3L']<-'3'
        startIndex<-which(qtldat$chr==peakChr & qtldat$Ppos==peakPos)
        BCI<-bayesint(qtldat,Bchr,startIndex,BCIprob)
        rownames(BCI)<-c('Lower Bound','Upper Bound')  
        allCI<-list(CIs,BCI)
        names(allCI)<-c('LODdrop','BCI')
        return(allCI)
        
      }else{
        stop("method must be either 'LODdrop', 'BCI', or 'both'")
      }
    }
  }
}

#the following function was modified from the qtl package
bayesint <-
  function(results, chr, qtl.index,prob)
  {    
    results <- results[results[,1]==chr & results[,'Gpos'] > results[qtl.index,'Gpos']-5 & results[,'Gpos'] < results[qtl.index,'Gpos']+5,]
    
    loc <- results[,2]
    width <- diff(( c(loc[1],loc) + c(loc, loc[length(loc)]) )/ 2)
    
    area <- 10^results[,'LOD']*width
    area <- area/sum(area)
    
    o <- order(results[,'LOD'], decreasing=TRUE)
    
    cs <- cumsum(area[o])
    
    wh <- min((1:length(loc))[cs >= prob])
    int <- range(o[1:wh])
    
    
    rn <- rownames(results)[c(int[1],o[1],int[2])]
    # look for duplicate rows
    if(any(table(rn)> 1)) {
      rn[2] <- paste(rn[2], "")
      if(rn[1] == rn[3]) rn[3] <- paste(rn[3], " ")
    }
    
    results <- results[c(int[1],int[2]),]
    return(results)
  }
