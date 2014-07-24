##' \code{entropy.pos} calculates the entropy (proportion of missing 
##' information) at a given position.
##' 
##' @title Entropy at a position
##'   
##' @aliases entropy.pos
##'   
##' @param peakChr character vector of length one. Must be one of the 
##'   major chromosome arms in the \emph{Drosophila} genome 
##'   ('X','2L','2R','3L',or '3R').
##'   
##' @param peakPos numeric vector of length one. A position in base 
##'   pairs in the DSPR position list (every 10kb).
##'   
##' @param design a character string. For inbred RIL designs:
##'   'inbredA', 'inbredB'. For cross designs: AAcross, BBcross, or
##'   'ABcross'. A and B refer to the pA and pB set of DSPR RILs.
##'   
##' @param phenotype.dat \code{data.frame} containing phenotype data. 
##'   For inbred designs, there must be a column of numeric RIL ids 
##'   (must be named patRIL). For the ABcross design, there must be
##'   both a patRIL and matRIL column specifying the pA and pB RIL
##'   ids. Cross designs also require a sex column for correct 
##'   specification of the genotypes on the X chromosome.
##'   
##' @param id.col a character string identifying the name of the 
##'   column containing unique ids for the samples. e.g. for an inbred
##'   design, the patRIL column can be used as the id.
##'   
##' @param sex a character string (either 'm' or 'f') specifying the 
##'   sex of the measured individuals. This argument must be supplied 
##'   for the AB cross design for correct specification of the
##'   genotypes on the X chromosome.
##'   
##' @return A numeric vector: the entropy at the given position for
##'   the set of RILs in the phenotype.dat \code{data.frame}. In the
##'   case of the ABcross, the A and B entropy are calculated
##'   separately and both are returned.
##'   
##' @author Elizabeth King (\email{kingeg@@missouri.edu})
##'   
##' @references Shannon, C.E. 1984. A mathematical theory of
##'   communication. \emph{Bell System Technical Journal} 27(3):
##'   379-423. 
##'   \url{http://en.wikipedia.org/wiki/Information_theory#Entropy}
##' 
##' @export
##' 
entropy.pos<-function(peakChr,peakPos,design,phenotype.dat,id.col,sex)
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
      data(list=objname,envir=environment())
    } else{
      con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                       objname, ".rda", sep = ""))
      load(con)
      close(con)
    }
    
    genotypes<-get(objname)
    rm(list=objname)
    
    genotypes<-merge(genotypes, phenotype.dat,by.x="ril",by.y="patRIL")
    mat.states<-as.matrix(genotypes[,c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                       "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                       "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                       "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
  }else{
    if(design=='inbredB')
    {
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
        data(list=objname,envir=environment())
      } else{
        con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataB/",
                         objname, ".rda", sep = ""))
        load(con)
        close(con)
      }
      
      genotypes<-get(objname)
      rm(list=objname)
      genotypes<-merge(genotypes, phenotype.dat,by.x="ril",by.y="patRIL")
      mat.states<-as.matrix(genotypes[,c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                         "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                         "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                         "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
    }else{
      if(design=='AAcross')
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
          data(list=objname,envir=environment())
        } else{
          con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                           objname, ".rda", sep = ""))
          load(con)
          close(con)
        }
        
        genotypes<-get(objname)
        rm(list=objname)
        
        patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
        matgeno<-merge(phenotype.dat,genotypes,by.x="matRIL",by.y="ril")
        patgeno<-patgeno[patgeno$id %in% matgeno$id,]
        matgeno<-matgeno[matgeno$id %in% patgeno$id,]
        patgeno<-patgeno[order(patgeno$id),]
        matgeno<-matgeno[order(matgeno$id),]
        
        
        
        if(peakChr=='X')
        {
          fpatgeno<-as.matrix(patgeno[patgeno$sex=='f',c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                                         "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                                         "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                                         "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
          fmatgeno<-as.matrix(matgeno[matgeno$sex=='f',c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                                         "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                                         "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                                         "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
          
          mmatgeno<-as.matrix(matgeno[matgeno$sex=='m',c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                                         "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                                         "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                                         "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
          
          mat.states<-rbind(fpatgeno,fmatgeno,mmatgeno)
          
        }else{
          patgenotypes<-as.matrix(patgeno[,c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                             "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                             "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                             "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
          matgenotypes<-as.matrix(matgeno[,c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                             "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                             "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                             "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])
          mat.states<-rbind(patgenotypes,matgenotypes)
        }
        
      }else{
        if(design=='BBcross')
        {
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
            data(list=objname,envir=environment())
          } else{
            con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataB/",
                             objname, ".rda", sep = ""))
            load(con)
            close(con)
          }
          genotypes<-get(objname)
          rm(list=objname)
          patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
          matgeno<-merge(phenotype.dat,genotypes,by.x="matRIL",by.y="ril")
          patgeno<-patgeno[patgeno$id %in% matgeno$id,]
          matgeno<-matgeno[matgeno$id %in% patgeno$id,]
          patgeno<-patgeno[order(patgeno$id),]
          matgeno<-matgeno[order(matgeno$id),]
          
          if(peakChr=='X')
          {
            fpatgeno<-as.matrix(patgeno[patgeno$sex=='f',c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                                           "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                                           "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                                           "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
            fmatgeno<-as.matrix(matgeno[matgeno$sex=='f',c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                                           "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                                           "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                                           "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
            
            mmatgeno<-as.matrix(matgeno[matgeno$sex=='m',c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                                           "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                                           "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                                           "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
            
            mat.states<-rbind(fpatgeno,fmatgeno,mmatgeno)
            
          }else{
            patgenotypes<-as.matrix(patgeno[,c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                               "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                               "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                               "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
            matgenotypes<-as.matrix(matgeno[,c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                               "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                               "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                               "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
            mat.states<-rbind(patgenotypes,matgenotypes)
          }
        }else{
          if(design=='ABcross')
          {
            if(require(DSPRqtlDataA)){
              
              use.package <- TRUE
            } else {
              message("Loading data from flyrils.org.\n
                      Consider installing DSPRqtlData[A/B] packages
                      for faster performance.\n")
              use.package <- FALSE
            }
            objnameA<-paste("A_",peakChr,"_",format(peakPos, sci = FALSE),sep="")
            
            if(use.package){
              data(list=objnameA,envir=environment())
            } else{
              con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                               objnameA, ".rda", sep = ""))
              load(con)
              close(con)
            }
            if(require(DSPRqtlDataB)){
              
              use.package <- TRUE
            } else {
              message("Loading data from flyrils.org.\n
                      Consider installing DSPRqtlData[A/B] packages
                      for faster performance.\n")
              use.package <- FALSE
            }
            objnameB<-paste("B_",peakChr,"_",format(peakPos, sci = FALSE),sep="")
            
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
            rm(list=objnameA)
            rm(list=objnameB)
            
            ABphenotype.dat<-phenotype.dat[phenotype.dat$patRIL<21000,]
            BAphenotype.dat<-phenotype.dat[phenotype.dat$patRIL>21000,]
            
            if(peakChr=='X' & sex=='m')
            {
              
              if(nrow(ABphenotype.dat)>0 & nrow(BAphenotype.dat)>0){stop("ABcross designs measuring males must all be a single type: 
                                                                         either A males to B females or B males to A females.")}
              
              if(nrow(BAphenotype.dat)>0)
              {
                matgeno<-merge(BAphenotype.dat,Agenotypes,by.x='matRIL',by.y='ril')
                matgeno<-merge(matgeno,Bgenotypes,by.x='patRIL',by.y='ril',sort=FALSE)
                mat.states<-as.matrix(matgeno[order(matgeno$id),c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                                                  "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                                                  "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                                                  "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8")])  
              }
              
              if(nrow(ABphenotype.dat)>0)
              {
                matgeno<-merge(ABphenotype.dat,Bgenotypes,by.x='matRIL',by.y='ril')
                matgeno<-merge(matgeno,Agenotypes,by.x='patRIL',by.y='ril') 
                mat.states<-as.matrix(matgeno[order(matgeno$id),c("B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                                                  "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                                                  "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                                                  "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
                
              }
              
              }else{           
                ABgenotypes<-merge(ABphenotype.dat,Agenotypes,by.x='patRIL',by.y='ril') 
                ABgenotypes<-merge(ABgenotypes, Bgenotypes, by.x='matRIL',by.y='ril',sort=FALSE)
                
                BAgenotypes<-merge(BAphenotype.dat,Agenotypes,by.x='matRIL',by.y='ril') 
                BAgenotypes<-merge(BAgenotypes, Bgenotypes, by.x='patRIL',by.y='ril',sort=FALSE)
                
                genotypes<-rbind(ABgenotypes,BAgenotypes)
                
                mat.states<-as.matrix(genotypes[,c("A1A1","A1A2","A1A3","A1A4","A1A5","A1A6","A1A7",
                                                   "A1A8","A2A2","A2A3","A2A4","A2A5","A2A6","A2A7","A2A8","A3A3","A3A4","A3A5",
                                                   "A3A6","A3A7","A3A8","A4A4","A4A5","A4A6","A4A7","A4A8","A5A5","A5A6","A5A7","A5A8",
                                                   "A6A6","A6A7","A6A8","A7A7","A7A8","A8A8",
                                                   "B1B1","B1B2","B1B3","B1B4","B1B5","B1B6","B1B7",
                                                   "B1B8","B2B2","B2B3","B2B4","B2B5","B2B6","B2B7","B2B8","B3B3","B3B4","B3B5",
                                                   "B3B6","B3B7","B3B8","B4B4","B4B5","B4B6","B4B7","B4B8","B5B5","B5B6","B5B7","B5B8",
                                                   "B6B6","B6B7","B6B8","B7B7","B7B8","B8B8")])
                
              }#else X/sex close 
          }else{
            stop("design must be one of: inbredA, inbredB, AAcross, BBcross, or ABcross")
          }#ABelse
            }#BBelse
        }#AAelse
      }#B else
    }#A else
  
  if(design=='ABcross' & ncol(mat.states)==72)
  {
    Asums<-apply(mat.states[,1:36],1,entropy)
    Bsums<-apply(mat.states[,37:72],1,entropy)
    Ainfo<-mean(Asums)
    Binfo<-mean(Bsums)
    info<-c(Ainfo,Binfo)
    names(info)<-c('A entropy','B entropy')
    return(info)
    
  }else{
    sums<-apply(mat.states,1,entropy)
    info<-mean(sums)
    return(info)
  }
  }


#undocumented functions for entropy.pos
entropy<-function(states)
  
{
  total<-0
  for (j in 1:length(states))
  {
    value<-(states[j]*eln(states[j]))/log(length(states))
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