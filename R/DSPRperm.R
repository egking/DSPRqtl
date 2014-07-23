##' \code{DSPRperm} performs a permutation test for a DSPR dataset.
##' 
##' The permutation test can take a very long time to run (over 24
##' hrs).
##' 
##' @title DSPR permutation test
##'   
##' @aliases DSPRperm print.pt
##'   
##' @param model an object of class formula: a symbolic description of
##'   the null model to be fitted at each position (e.g., 
##'   \code{phenotype ~ 1}). The genotype effects to be fitted will 
##'   be added based on \code{design}.
##'   
##' @param design a character string. One of either 'inbredA', 
##'   'inbredB', or 'ABcross' corresponding to the pA and pB set of 
##'   inbred RILs or the pA-pB cross design. For round robin designs 
##'   or other cross designs, use the more flexible DSPRgenos and 
##'   standard model fitting functions in R.
##'   
##' @param phenotype.dat \code{data.frame} containing phenotype data. 
##'   For inbred designs, there must be a column of numeric RIL ids 
##'   (must be named patRIL). For the ABcross design, there must be
##'   both a patRIL and matRIL column specifying the pA and pB RIL
##'   ids.
##'   
##' @param id.col a character string identifying the name of the 
##'   column containing unique ids for the samples. e.g. for an inbred
##'   design, the patRIL column can be used as the id.
##'   
##' @param batch A numeric vector of length one specifying the number 
##'   of positions to be examined at a time. A larger number will use 
##'   more memory but can be faster. Default is 1000.
##'   
##' @param niter A numeric vector of length one specifying the number 
##'   of permutations to run. Default is 1000.
##'   
##' @param alpha The alpha level for the genome-wide signficance 
##'   threshold. Default is 0.05. Raw maximum LOD scores are also 
##'   returned and \code{\link{quantile}} can be used to test other 
##'   alphas.
##'   
##' @param sex a character string (either 'm' or 'f') specifying the 
##'   sex of the measured individuals. This argument must be supplied 
##'   for a cross design for correct specification of the genotypes on
##'   the X chromosome.
##' 
##' @return A list of class \code{pt} containing:
##' \item{maxLODs}{A vector containing the maximum LOD score obtained
##' for each permutation of the data.}
##' \item{alpha}{the specified alpha level}
##' \item{threshold}{the signficance threshold at the specified alpha
##' level}
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##'
DSPRperm<-function(model,design,phenotype.dat,id.col, batch=1000,niter=1000,alpha=0.05,sex)
{ 
  #CHECK THAT ABCROSS HAS SPECIFIED SEX
  if(missing(sex)){
    if(design=='ABcross')
    {
      stop("If using the ABcross design, you must specify the 
           sex of the offspring so the genotypes on the X chromosome 
           can be specified correctly")
    }else{
      sex<-NA
    }
    
  }
  
  ##ASSIGN ID COLUMN
  phenotype.dat$id<-phenotype.dat[,id.col]
  
  ##CHECK FOR REQUIRED COLUMNS
  if(design=='inbredA'|design=='inbredB')
  {
    
    if(!('patRIL' %in% colnames(phenotype.dat)))
    {
      stop("phenotype data frame must contain a patRIL column")
    }
  }
  
  if(design=='ABcross'|design=='AAcross'|design=='BBcross')
  {
    if(!('patRIL' %in% colnames(phenotype.dat)) | 
         !('matRIL' %in% colnames(phenotype.dat)) |
         !('sex' %in% colnames(phenotype.dat)))
    {
      stop("for cross designs, phenotype data frame must 
           contain the following columns: patRIL, matRIL, and sex")
    }
  }
  
  
  
  if(design=='inbredA'|design=='inbredB'|design=='ABcross')
  {
    if(design=='inbredA'|design=='ABcross')
    {
      if(require(DSPRqtlDataA)){
        
        use.package <- TRUE
      } else {
        message("Loading data from flyrils.org.\n
                Consider installing DSPRqtlData[A/B] packages
                for faster performance.\n")
        use.package <- FALSE
      }
    }
    
    if(design=='inbredB'|design=='ABcross')
    {
      if(require(DSPRqtlDataB)){
        
        use.package <- TRUE
      } else {
        message("Loading data from flyrils.org.\n
                Consider installing DSPRqtlData[A/B] packages
                for faster performance.\n")
        use.package <- FALSE
      }
    }
  }else{
    stop("must specify valid design: inbredA, inbredB, or ABcross")
  }
    
  #get list of positions
  data(positionlist_wgenetic,envir=environment())
  
  
  
   
  if(design=='inbredA'|design=='inbredB')
  {
    
    if(design=='inbredA')
    {
      objname<-paste("A_",poslist[1,1],"_",format(poslist[1,2], sci = FALSE),sep="")
      
    }else{
      objname<-paste("B_",poslist[1,1],"_",format(poslist[1,2], sci = FALSE),sep="")
    }
    data(list=objname,envir=environment())
    genotypes<-get(objname)
    #this makes sure all your phenotyped rils are in the matrix of genotypes
    phenotype.dat<-phenotype.dat[phenotype.dat$patRIL %in% genotypes$ril,]
    #order phenotype.dat by id
    phenotype.dat<-phenotype.dat[order(phenotype.dat$id),]
    rm(list=objname)
    
  }else{
    objnameA<-paste("A_",poslist[1,1],"_",format(poslist[1,2], sci = FALSE),sep="")
    objnameB<-paste("B_",poslist[1,1],"_",format(poslist[1,2], sci = FALSE),sep="")
    data(list=objnameA,envir=environment())
    Agenotypes<-get(objnameA)
    data(list=objnameB,envir=environment())
    Bgenotypes<-get(objnameB)
    ABphenotype.dat<-phenotype.dat[phenotype.dat$patRIL<21000,]
    BAphenotype.dat<-phenotype.dat[phenotype.dat$patRIL>21000,]
    ABphenotype.dat<-ABphenotype.dat[ABphenotype.dat$patRIL %in% Agenotypes$ril & ABphenotype.dat$matRIL %in% Bgenotypes$ril,]
    BAphenotype.dat<-BAphenotype.dat[BAphenotype.dat$patRIL %in% Bgenotypes$ril & BAphenotype.dat$matRIL %in% Agenotypes$ril,]
    phenotype.dat<-rbind(ABphenotype.dat,BAphenotype.dat)
    #order phenotype.dat by id
    phenotype.dat<-phenotype.dat[order(phenotype.dat$id),]
    rm(list=objnameA)
    rm(list=objnameB)    
  }
  
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
            data(list=objname,envir=environment())
          } else{
            con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                             objname, ".rda", sep = ""))
            load(con)
            close(con)
          }
          genotypes<-get(objname)
          patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
          genos<-as.matrix(patgeno[order(patgeno$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
          row.names(genos)<-patgeno$id
          big.list[[counter]]<-genos
          rm(list=objname)
          counter<-counter+1
          #time2<-Sys.time()
        }else{
          objname<-paste("B_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
          
          if(use.package){
            data(list=objname,envir=environment())
          } else{
            con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataB/",
                             objname, ".rda", sep = ""))
            load(con)
            close(con)
          }
          
          genotypes<-get(objname)
          patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
          genos<-as.matrix(patgeno[order(patgeno$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
          row.names(genos)<-patgeno$id
          big.list[[counter]]<-genos
          rm(list=objname)
          counter<-counter+1
        }#B else close
        
      }else{
        objnameA<-paste("A_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
        
        if(use.package){
          data(list=objnameA,envir=environment())
        } else{
          con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                           objnameA, ".rda", sep = ""))
          load(con)
          close(con)
        }
        objnameB<-paste("B_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
        
        if(use.package){
          data(list=objnameB,envir=environment())
        } else{
          con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataB/",
                           objnameB, ".rda", sep = ""))
          load(con)
          close(con)
        }
        Agenotypes<-get(objnameA)
        Bgenotypes<-get(objnameB)
        
        ABphenotype.dat<-phenotype.dat[phenotype.dat$patRIL<21000,]
        BAphenotype.dat<-phenotype.dat[phenotype.dat$patRIL>21000,]
        
        if(poslist[i,1]=='X' & sex=='m')
        {
          
          if(nrow(ABphenotype.dat)>0 & nrow(BAphenotype.dat)>0){stop("ABcross designs measuring males must all be a single type: 
                                                                     either A males to B females or B males to A females.")}
          
          if(nrow(BAphenotype.dat)>0)
          {
            matgeno<-merge(BAphenotype.dat,Agenotypes,by.x='matRIL',by.y='ril')
            matgeno<-merge(matgeno,Bgenotypes,by.x='patRIL',by.y='ril',sort=FALSE)
            genos<-as.matrix(matgeno[order(matgeno$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
            row.names(genos)<-matgeno$id
            big.list[[counter]]<-genos
            rm(list=objnameA)
            rm(list=objnameB)
            counter<-counter+1
            
          }
          
          if(nrow(ABphenotype.dat)>0)
          {
            matgeno<-merge(ABphenotype.dat,Bgenotypes,by.x='matRIL',by.y='ril') 
            matgeno<-merge(matgeno,Agenotypes,by.x='patRIL',by.y='ril') 
            genos<-as.matrix(matgeno[order(matgeno$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
            row.names(genos)<-matgeno$id
            big.list[[counter]]<-genos
            rm(list=objnameA)
            rm(list=objnameB)
            counter<-counter+1
          }
          
          }else{           
            ABgenotypes<-merge(ABphenotype.dat,Agenotypes,by.x='patRIL',by.y='ril') 
            ABgenotypes<-merge(ABgenotypes, Bgenotypes, by.x='matRIL',by.y='ril',sort=FALSE)
            
            BAgenotypes<-merge(BAphenotype.dat,Agenotypes,by.x='matRIL',by.y='ril') 
            BAgenotypes<-merge(BAgenotypes, Bgenotypes, by.x='patRIL',by.y='ril',sort=FALSE)
            
            genotypes<-rbind(ABgenotypes,BAgenotypes)
            genotypes<-genotypes[order(genotypes$id),]
            
            genos<-as.matrix(genotypes[,c("AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8",
                                          "BB1","BB2","BB3","BB4","BB5","BB6","BB7","BB8")])
            row.names(genos)<-genotypes$id 
            big.list[[counter]]<-genos
            rm(list=objnameA)
            rm(list=objnameB)
            counter<-counter+1
            
          }#else X/sex close
        
      }# AB else close
      
    }#i close
    
    
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
    full.lod.set[pstart:pend,]<-all.lms-L.n
   }else{
    
     for(jj in 1:niter)
    {
      pheno<-phenotype.dat[ysamps[,jj],]
      
    
    #get model likelihoods at each position (will take several minutes)
    lms<-lapply(big.list, function(x) LL.alt(x,model,pheno)) 
      
   
      
    #get LOD scores
    lod.set<-unlist(lms)-L.n
    
    full.lod.set[pstart:pend,jj]<-lod.set  
  
    }#niter close
   }#else close
    rm(big.list)
    
  }#batch close
  rm(poslist)

 if(design=='ABcross' & sex=='m')
 {
   maxlodsX<-apply(full.lod.set[poslist$chr=='X',],2,max)
   maxlodsA<-apply(full.lod.set[poslist$chr!='X',],2,max)
   mm<-list(maxlodsX,maxlodsA)
   names(mm)<-c('X','Autosomes')
   th<-list(quantile(maxlodsX,1-alpha),quantile(maxlodsA,1-alpha))
   names(th)<-c('X','Autosomes')
   perm.list<-list("maxLODs"=mm,
                   "alpha"=alpha,
                   "threshold"=th
   )
   class(perm.list)<-'pt'
   return(perm.list)
 }else{
  maxlods<-apply(full.lod.set,2,max)
  
  perm.list<-list("maxLODs"=maxlods,
                  "alpha"=alpha,
                  "threshold"=quantile(maxlods,1-alpha)
                  )
  class(perm.list)<-'pt'
  return(perm.list)
}      
} #function close

print.pt <- function(x, digits = 4, ...){
  cat("\n\n")
  cat("Genome-wide signficance threshold (alpha = ",x$alpha,"):",x$threshold,"\n\n\n")
  }
  
