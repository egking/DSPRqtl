##' \code{DSPRperm} performs a permutation test for a DSPR dataset.
##'
##' The permutation test can take a very long time to run (over 24 hrs).
##' 
##' @title DSPR permutation test
##' 
##' @aliases DSPRperm print.pt
##' 
##' 
##' @param model an object of class formula: a symbolic description of the 
##' null model to be fitted at each position (e.g., \code{phenotype \~ 1}).  
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
##' @param niter A numeric vector of length one specifying the number of permutations 
##' to run. Default is 1000.
##' 
##' @param alpha The alpha level for the genome-wide signficance threshold. Default is 0.05. 
##' Raw maximum LOD scores are also returned and \code{\link{quantile}} can be used to test 
##' other alphas.
##' 
##' @return A list of class \code{pt} containing:
##' \item{maxLODs}{A vector containing the maximum LOD score obtained for each permutation of 
##' the data.} 
##' \item{alpha}{the specified alpha level}
##' \item{threshold}{the signficance threshold at the specified alpha level}
##' 
##' 
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##'
##' @S3method print pt
DSPRperm<-function(model,design,phenotype.dat, batch=1000,niter=1000,alpha=0.05)
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
  
  #order phenotype.dat by ril
  phenotype.dat<-phenotype.dat[order(phenotype.dat$patRIL),]
  
  if(design=='inbredA'|design=='inbredB')
  {
    
    if(design=='inbredA')
    {
      objname<-paste("A_",poslist[1,1],"_",format(poslist[1,2], sci = FALSE),sep="")
      
    }else{
      objname<-paste("B_",poslist[1,1],"_",format(poslist[1,2], sci = FALSE),sep="")
    }
  }
  data(list=objname)
  genotypes<-get(objname)
  #this makes sure all your phenotyped rils are in the matrix of genotypes
  phenotype.dat<-phenotype.dat[phenotype.dat$patRIL %in% genotypes$ril,]
  rm(list=objname,pos=.GlobalEnv)
  
  #set up index to sample
  index<-seq(1:nrow(phenotype.dat))
  ysamps<-matrix(,nrow(phenotype.dat),niter)
  for (j in 1:niter)
  {
    ysamps[,j]<-sample(index)
  }
  
  if(batch=='full'){batch=nrow(poslist)}
  
  batches<-seq(1,nrow(poslist),by=batch)
  
  full.lod.set<-matrix(,nrow(poslist),niter)
  for (b in 1:length(batches))
  {
    start<-batches[b]  
    if(b==length(batches)){end<-nrow(poslist)}else{end<-batches[b+1]-1} 
    
    big.list<-vector('list',(end-start)+1)
    counter<-1
    
    #prepare data
    for (i in start:end) 
    {  
      
      if(design=='inbredA'|design=='inbredB')
      {
        
        if(design=='inbredA')
        {
          objname<-paste("A_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
          data(list=objname)
          genotypes<-get(objname)
          patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
          genos<-as.matrix(patgeno[order(patgeno$patRIL),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
          row.names(genos)<-patgeno$patRIL
          big.list[[counter]]<-genos
          counter<-counter+1
          rm(list=objname,pos=.GlobalEnv)
          
        }else{
          
          objname<-paste("B_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
          data(list=objname)
          genotypes<-get(objname)
          patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
          genos<-as.matrix(patgeno[order(patgeno$patRIL),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
          row.names(genos)<-patgeno$patRIL
          big.list[[counter]]<-genos
          counter<-counter+1
          rm(list=objname,pos=.GlobalEnv)
          
        }
        
      }
      
    }
    
    
    #get null likelihood
    null.mod<-lm(model,data=phenotype.dat)
    L.n<-logLik(null.mod)/log(10)
    
    
    if(grepl("~\\s*1\\s*$", deparse(model))){
    pheno<-matrix(,nrow(phenotype.dat),niter)
    for (jj in 1:niter)
    {
      
      pheno[,jj]<-phenotype.dat[ysamps[,jj],as.character(model)[2]]
    }
    
    lms<-lapply(big.list, function(x) LL.multi(x,model,pheno))
    
    
    all.lms <- array(dim=c(length(lms),dim(lms[[1]])[1]))
    
    for(zz in seq(along=lms)) all.lms[zz,] <- lms[[zz]]
    full.lod.set[start:end,]<-all.lms-L.n
   }else{
    
     for(jj in 1:niter)
    {
      pheno<-phenotype.dat[ysamps[,jj],]
      
    
    #get model likelihoods at each position (will take several minutes)
    lms<-lapply(big.list, function(x) LL.alt(x,model,pheno)) 
      
   
      
    #get LOD scores
    lod.set<-unlist(lms)-L.n
    
    full.lod.set[start:end,j]<-lod.set  
  
    }#niter close
   }#else close
    rm(big.list)
    
  }#batch close
  rm(poslist,pos=.GlobalEnv)
  
  maxlods<-apply(full.lod.set,2,max)
  
  perm.list<-list("maxLODs"=maxlods,
                  "alpha"=alpha,
                  "threshold"=quantile(maxlods,1-alpha)
                  )
  class(perm.list)<-'pt'
  return(perm.list)
        
} #function close

print.pt <- function(x, digits = 4, ...){
  cat("\n\n")
  cat("Genome-wide signficance threshold (alpha = ",x$alpha,"):",x$threshold,"\n\n\n")
  }
  
