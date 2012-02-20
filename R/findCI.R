findCI<-function(peakChr,peakPos,qtldat,peak,drop)
{
  startIndex<-which(qtldat$chr==peakChr & qtldat$Ppos==peakPos)
  l.ci<-peak-drop  #drop threshold
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
  return()
}