##' Function to perform a genome scan for data generated with the 
##' DSPR RILs.
##'
##' 
##' 
##' @title DSPR Genome Scan
##' 
##' @param model an object of class formula: a symbolic description of the 
##' null model to be fitted at each position (e.g., \code{phenotype ~ 1}).  
##' The genotype effects to be fitted 
##' will be added based on \code{design}. 
##' 
##' @param design a character string. One of either 'inbredA' or 'inbredB'
##' corresponding to the pA and pB set of inbred RILs. Other crossing designs 
##' will be supported in the future. 
##' 
##' @param phenotype.dat \code{data.frame} containing a column of 
##' ril ids (must be named patRIL) and phenotypes.
##' 
##' @param batch A numeric vector of length one specifying the number of positions 
##' to be examined at a time. A larger number will use more memory but can be faster. 
##' Default is 1000.
##' 
##' @return A list of class \code{gscan} containing:
##' \item{LODscores}{A \code{data.frame} consisting of the chromosome, physical position (bp),
##' genetic position (cM) and LOD score for each significant peak.}
##' \item{model}{the model specification}
##' \item{design}{the design specification}
##' \item{phenotype}{the phenotype \code{data.frame} specified}
##' 
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##'
DSPRscan<-function(model,design,phenotype.dat, batch=1000)
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
    }else{
      if(require(DSPRqtlDataB)){
        
        use.package <- TRUE
      } else {
        message("Loading data from flyrils.org.\n
                Consider installing DSPRqtlData[A/B] packages
                for faster performance.\n")
        use.package <- FALSE
      }
    }
    
  }
  
  
  #get list of positions
  data(positionlist_wgenetic)
  #poslist<-poslist[1:5,]
  if(batch=='full'){batch=nrow(poslist)}
  
  batches<-seq(1,nrow(poslist),by=batch)
  
  full.lod.set<-numeric(length=nrow(poslist))
  for (b in 1:length(batches))
  {
    pstart<-batches[b]  
    if(b==length(batches)){pend<-nrow(poslist)}else{pend<-batches[b+1]-1} 
    
    big.list<-vector('list',(pend-pstart)+1)
    counter<-1
    #prepare data
  for (i in pstart:pend) 
  {  
    
    if(design=='inbredA'|design=='inbredB')
    {
      
      if(design=='inbredA')
      {
       # time1<-Sys.time()
        objname<-paste("A_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
        if(use.package){
          data(list=objname)
        } else{
          con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                           objname, ".rda", sep = ""))
          load(con)
          close(con)
          }
        genotypes<-get(objname)
        patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
        genos<-as.matrix(patgeno[order(patgeno$patRIL),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
        row.names(genos)<-patgeno$patRIL
        big.list[[counter]]<-genos
        rm(list=objname,pos=.GlobalEnv)
        counter<-counter+1
        #time2<-Sys.time()
        }else{
          objname<-paste("B_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
          
          if(use.package){
            data(list=objname)
          } else{
            con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataB/",
                             objname, ".rda", sep = ""))
            load(con)
            close(con)
          }

          genotypes<-get(objname)
          patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
          genos<-as.matrix(patgeno[order(patgeno$patRIL),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
          row.names(genos)<-patgeno$patRIL
          big.list[[counter]]<-genos
          rm(list=objname,pos=.GlobalEnv)
          counter<-counter+1
      }
      
    }
    
    
    
  }
  
  #order phenotype.dat by ril
  phenotype.dat<-phenotype.dat[order(phenotype.dat$patRIL),]
  #this makes sure all your phenotyped rils are in the matrix of genotypes
  phenotype.dat<-phenotype.dat[phenotype.dat$patRIL %in% rownames(big.list[[1]]),]
  
  #get null likelihood
  null.mod<-lm(model,data=phenotype.dat)
  L.n<-logLik(null.mod)/log(10)
  
  #get model likelihoods at each position (will take several minutes)
  lms<-lapply(big.list, function(x) LL.alt(x,model,phenotype.dat)) 
  rm(big.list)
  #get LOD scores
  lod.set<-unlist(lms)-L.n
  
  full.lod.set[pstart:pend]<-lod.set  
  }#batch close
  
  #put position and LOD scores in a data frame
  #qtl.results<-data.frame('chr'=poslist$chr,'Ppos'=poslist$Ppos,'Gpos'=poslist$Gpos,'LOD'=lod.set)
  qtl.results<-list("LODscores"=data.frame('chr'=poslist$chr,'Ppos'=poslist$pos,'Gpos'=poslist$Gpos,'LOD'=full.lod.set),
                    "model"=model,
                    "design"=design,
                    "phenotype"=phenotype.dat
                    )
 
  class(qtl.results)<-'gscan'
  return(qtl.results)
} #function close

