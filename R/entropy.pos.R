##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param peakChr 
##' @param peakPos 
##' @param phenotype.dat 
##' @return 
##' @author Elizabeth King
entropy.pos<-function(peakChr,peakPos,phenotype.dat)
{
  if(design=='inbredA'|design=='inbredB')
  {
    
    if(design=='inbredA')
    {
      if(require(DSPRqtlDataA)){
        
        use.package <- TRUE
      } else {
        message("Loading data from flyrils.org.\n
                Consider installing DSPRqtlData[A/B] packages
                for faster performance.\n")
        use.package <- FALSE
      }
      objname<-paste("A_",peakChr,"_",format(peakPos, sci = FALSE),sep="")
      
      if(use.package){
        data(list=objname)
      } else{
        con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                         objname, ".rda", sep = ""))
        load(con)
        close(con)
      }
      
      genotypes<-get(objname)
      rm(list=objname,pos=.GlobalEnv)
      
      genotypes<-merge(genotypes, phenotype.dat,by.x="ril",by.y="patRIL")
      mat.states<-as.matrix(genotypes[,c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                         "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                         "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                         "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
    }else{
      if(require(DSPRqtlDataB)){
        
        use.package <- TRUE
      } else {
        message("Loading data from flyrils.org.\n
                Consider installing DSPRqtlData[A/B] packages
                for faster performance.\n")
        use.package <- FALSE
      }
      
      objname<-paste("B_",peakChr,"_",format(peakPos, sci = FALSE),sep="")
      
      if(use.package){
        data(list=objname)
      } else{
        con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataB/",
                         objname, ".rda", sep = ""))
        load(con)
        close(con)
      }
      
      genotypes<-get(objname)
      rm(list=objname,pos=.GlobalEnv)
      genotypes<-merge(genotypes, phenotype.dat,by.x="ril",by.y="patRIL")
      mat.states<-as.matrix(genotypes[,c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                         "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                         "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                         "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
    }
  }
  
  
  sums<-apply(mat.states,1,entropy)
  info<-mean(sums)
  return(info)
  
}


#undocumented functions for entropy.pos
entropy<-function(states)
  
{
  total<-0
  for (j in 1:length(states))
  {
    value<-(states[j]*eln(states[j]))/log(8)
    total<-total+value  
  }
  return(-total)
}

eln<-function(log_value)
{
  if (log_value==0)
  {
    return(0)	
  }else{
    log_value<-log(log_value)
    return(log_value)
  }
  
}
