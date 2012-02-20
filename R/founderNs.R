##' Numbers of founder genotypes
##'
##' \code{founderNs} calculates the number of founder genotypes represented in the DSPR
##' RILs at a given position.
##' 
##' @title Numbers of founder genotypes
##' @param peakChr character vector of length one. 
##' Must be one of the major chromosome arms in the 
##' \emph{Drosophila} genome ('X','2L','2R','3L',or '3R').
##' 
##' @param peakPos numeric vector of length one. A position in 
##' base pairs in the DSPR position list (every 10kb). See 
##' \code{\link{positionlist_wgenetic}}.
##' 
##' @param design a character string. One of either 'inbredA' or 'inbredB'
##' corresponding to the pA and pB set of inbred RILs. Other crossing designs 
##' will be supported in the future. 
##' 
##' @param phenotype.dat \code{data.frame} containing a column of 
##' ril ids (must be named patRIL) and phenotypes.
##'
##' @return A named numeric vector consisting of the numbers of RILs in the 
##' phenotype.dat \code{data.frame} with each founder genotype at the given position. 
##' A RIl is assigned a founder genotype if the HMM probability is over 0.95. The number 
##' of RILs with a heterozygous genotype and an uncertain genotype are also returned. 
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##' 

founderNs<-function(peakChr,peakPos,design,phenotype.dat)
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
      allgenotypes<-as.matrix(genotypes[,c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                           "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                           "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                           "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
      max.probs<-apply(allgenotypes,1,max)
      n.under<-length(max.probs[max.probs<0.95])
      A1.n<-nrow(allgenotypes[allgenotypes[,'A1A1']>=0.95,])
      A2.n<-nrow(allgenotypes[allgenotypes[,'A2A2']>=0.95,])
      A3.n<-nrow(allgenotypes[allgenotypes[,'A3A3']>=0.95,])
      A4.n<-nrow(allgenotypes[allgenotypes[,'A4A4']>=0.95,])
      A5.n<-nrow(allgenotypes[allgenotypes[,'A5A5']>=0.95,])
      A6.n<-nrow(allgenotypes[allgenotypes[,'A6A6']>=0.95,])
      A7.n<-nrow(allgenotypes[allgenotypes[,'A7A7']>=0.95,])
      A8.n<-nrow(allgenotypes[allgenotypes[,'A8A8']>=0.95,])
      max.het.probs<-apply(allgenotypes[,c('A1A2','A1A3','A1A4','A1A5','A1A6','A1A7','A1A8',
                                           'A2A3','A2A4','A2A5','A2A6','A2A7','A2A8',
                                           'A3A4','A3A5','A3A6','A3A7','A3A8',
                                           'A4A5','A4A6','A4A7','A4A8',
                                           'A5A6','A5A7','A5A8',
                                           'A6A7','A6A8',
                                           'A7A8'
                                           )],1,max)
      Ah.n<-length(max.het.probs[max.het.probs>=0.95])
      Ns<-c(A1.n,A2.n,A3.n,A4.n,A5.n,A6.n,A7.n,A8.n,Ah.n,n.under)
      names(Ns)<-c('A1','A2','A3','A4','A5','A6','A7','A8','Hets','Uncertain')
      
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
      allgenotypes<-as.matrix(genotypes[,c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                           "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                           "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                           "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
      max.probs<-apply(allgenotypes,1,max)
      n.under<-length(max.probs[max.probs<0.95])
      B1.n<-nrow(allgenotypes[allgenotypes[,'B1B1']>=0.95,])
      B2.n<-nrow(allgenotypes[allgenotypes[,'B2B2']>=0.95,])
      B3.n<-nrow(allgenotypes[allgenotypes[,'B3B3']>=0.95,])
      B4.n<-nrow(allgenotypes[allgenotypes[,'B4B4']>=0.95,])
      B5.n<-nrow(allgenotypes[allgenotypes[,'B5B5']>=0.95,])
      B6.n<-nrow(allgenotypes[allgenotypes[,'B6B6']>=0.95,])
      B7.n<-nrow(allgenotypes[allgenotypes[,'B7B7']>=0.95,])
      B8.n<-nrow(allgenotypes[allgenotypes[,'B8B8']>=0.95,])
      max.het.probs<-apply(allgenotypes[,c('B1B2','B1B3','B1B4','B1B5','B1B6','B1B7','B1B8',
                                           'B2B3','B2B4','B2B5','B2B6','B2B7','B2B8',
                                           'B3B4','B3B5','B3B6','B3B7','B3B8',
                                           'B4B5','B4B6','B4B7','B4B8',
                                           'B5B6','B5B7','B5B8',
                                           'B6B7','B6B8',
                                           'B7B8'
                                           )],1,max)
      Bh.n<-length(max.het.probs[max.het.probs>=0.95])
      Ns<-c(B1.n,B2.n,B3.n,B4.n,B5.n,B6.n,B7.n,B8.n,Bh.n,n.under)
      names(Ns)<-c('B1','B2','B3','B4','B5','B6','B7','B8','Hets','Uncertain')
      
    }
  }
  return(Ns)
}
