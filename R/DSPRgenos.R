##' Function to generate genotype probabilities across the genome for 
##' a given DSPR dataset. This function is most useful for those wishing 
##' to fit their own models. 
##'
##' @title DSPR genotype probabilities
##' 
##' @param design a character string. For inbred RIL designs: 'inbredA', 
##'   'inbredB'. For cross designs: AAcross, BBcross, or 'ABcross'. 
##'   A and B refer to the pA and pB set of 
##'   DSPR RILs.   
##'   
##' @param phenotype.dat \code{data.frame} containing phenotype data. 
##' For inbred designs, there must be a column of numeric RIL ids 
##' (must be named patRIL). For cross designs, there must be both
##' a patRIL and matRIL column specifying the maternal and paternal 
##' RIL ids. Cross designs also require a sex column for correct 
##' specification of the genotypes on the X chromosome. 
##'
##' @param id.col a character string identifying the name of the 
##' column containing unique ids for the samples. e.g. for an inbred
##' design, the patRIL column can be used as the id. 
##' 
##' @return A list containing:
##' \item{genolist}{a list containing the matrix of additive genotype 
##' probabilities at each position in the genome. The list is in the 
##' same order as the list of positions described below. Column names 
##' are the different DSPR haplotypes and row names are the unique ids 
##' provided in id.col}
##' \item{positions}{a \code{data.frame} containing regularly spaced 
##' positions (every 10KB) in the genome where the genotype probabilities 
##' are calculated. Columns are: chr = chromosome, Ppos = physical position 
##' (zero offset), Gpos = genetic position, Gaxis = cummulative genetic 
##' position.}
##' \item{phenotype}{the phenotype.dat \code{data.frame} ordered by the 
##' specified id column and in the same order as the genotype information
##' at each position in the genome in the geno list described above}
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##'


DSPRgenos<-function(design,phenotype.dat,id.col)
{
  
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
  
  
  ## CHECK IF DATA PACKAGES ARE INSTALLED AND LOAD THEM
  if(design=='inbredA'|design=='ABcross'|design=='AAcross')
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
    if(design=='inbredB'|design=='ABcross'|design=='BBcross')
    {
      if(require(DSPRqtlDataB)){
        
        use.package <- TRUE
        
      } else {
        message("Loading data from flyrils.org.\n
                Consider installing DSPRqtlData[A/B] packages
                for faster performance.\n")
        use.package <- FALSE
      }
    }else{
      stop("Design not valid. 
           Must be one of inbredA, inbredB, ABcross, 
           AAcross, or BBcross")
    }
    
    }
  
  #GET POSITION DATA 
  data(positionlist_wgenetic)
  
  #SET UP LIST FOR GENOTYPE MATRICIES
  big.list<-vector('list',nrow(poslist))
  
  if(design=='inbredA')
  {
    
    for (i in 1:length(big.list)) 
    {  
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
      genos<-as.matrix(patgeno[order(patgeno$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
      row.names(genos)<-patgeno$id
      big.list[[i]]<-genos
      rm(list=objname,pos=.GlobalEnv)
    }
    
  }else{
    if(design=='inbredB')
    {
      for (i in 1:length(big.list)) 
      {
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
        genos<-as.matrix(patgeno[order(patgeno$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
        row.names(genos)<-patgeno$id
        big.list[[i]]<-genos
        rm(list=objname,pos=.GlobalEnv)
        
      }  
    }else{
      if(design=='AAcross')
      {
        
        for (i in 1:length(big.list)) 
        {  
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
          matgeno<-merge(phenotype.dat,genotypes,by.x="matRIL",by.y="ril")
          patgeno<-patgeno[patgeno$id %in% matgeno$id,]
          matgeno<-matgeno[matgeno$id %in% patgeno$id,]
          patgeno<-patgeno[order(patgeno$id),]
          matgeno<-matgeno[order(matgeno$id),]
          
          if(poslist[i,1]=='X')
          {
            fpatgeno<-patgeno[patgeno$sex=='f',]
            fmatgeno<-matgeno[matgeno$sex=='f',]
            fgeno<-(fpatgeno[,c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')]
                    +fmatgeno[,c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])/2
            
            fgeno$id<-fpatgeno$id
            
            mmatgeno<-matgeno[matgeno$sex=='m',]
            mgeno<-mmatgeno[,c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8','id')]
            genos<-rbind(fgeno,mgeno)
            genos<-as.matrix(genos[order(genos$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
            
          }else{
            genos<-(as.matrix(patgeno[order(patgeno$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
                    +as.matrix(matgeno[order(matgeno$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')]))/2
          }
          
          row.names(genos)<-patgeno$id
          big.list[[i]]<-genos
          rm(list=objname,pos=.GlobalEnv)
        }
        
      }else{
        if(design=='BBcross')
        {
          for (i in 1:length(big.list)) 
          {
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
            matgeno<-merge(phenotype.dat,genotypes,by.x="matRIL",by.y="ril")
            patgeno<-patgeno[patgeno$id %in% matgeno$id,]
            matgeno<-matgeno[matgeno$id %in% patgeno$id,]
            patgeno<-patgeno[order(patgeno$id),]
            matgeno<-matgeno[order(matgeno$id),]
            
            if(poslist[i,1]=='X')
            {
              fpatgeno<-patgeno[patgeno$sex=='f',]
              fmatgeno<-matgeno[matgeno$sex=='f',]
              fgeno<-(fpatgeno[,c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')]
                      +fmatgeno[,c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])/2
              
              fgeno$id<-fpatgeno$id
              mmatgeno<-matgeno[matgeno$sex=='m',]
              mgeno<-mmatgeno[,c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8','id')]
              genos<-rbind(fgeno,mgeno)
              genos<-as.matrix(genos[order(genos$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
              
            }else{
              genos<-(as.matrix(patgeno[order(patgeno$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
                      +as.matrix(matgeno[order(matgeno$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')]))/2
            }
            
            row.names(genos)<-patgeno$id
            big.list[[i]]<-genos
            rm(list=objname,pos=.GlobalEnv)          
          }
        }else{
          
          if(design=='ABcross')
          {
            for (i in 1:length(big.list)) 
            {
              objnameA<-paste("A_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
              
              if(use.package){
                data(list=objnameA)
              } else{
                con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataA/",
                                 objnameA, ".rda", sep = ""))
                load(con)
                close(con)
              }
              objnameB<-paste("B_",poslist[i,1],"_",format(poslist[i,2], sci = FALSE),sep="")
              
              if(use.package){
                data(list=objnameB)
              } else{
                con <- url(paste("http://wfitch.bio.uci.edu/R/DSPRqtlDataB/",
                                 objnameB, ".rda", sep = ""))
                load(con)
                close(con)
              }
              Agenotypes<-get(objnameA)
              Bgenotypes<-get(objnameB)
              
              ABphenotype.dat<-phenotype.dat[phenotype.dat$patRIL<21000,]
              BAphenotype.dat<-phenotype.dat[phenotype.dat$patRI>21000,]
              
              ABgenotypes<-merge(ABphenotype.dat,Agenotypes,by.x='patRIL',by.y='ril') 
              ABgenotypes<-merge(ABgenotypes, Bgenotypes, by.x='matRIL',by.y='ril',sort=FALSE)
              
              BAgenotypes<-merge(BAphenotype.dat,Agenotypes,by.x='matRIL',by.y='ril') 
              BAgenotypes<-merge(BAgenotypes, Bgenotypes, by.x='patRIL',by.y='ril',sort=FALSE)
              
              genotypes<-rbind(ABgenotypes,BAgenotypes)
              genotypes<-genotypes[order(genotypes$id),]
              genos<-as.matrix(genotypes[,c("AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8","BB1","BB2","BB3","BB4","BB5","BB6","BB7","BB8")])
              row.names(genos)<-genotypes$id 
              
              if(poslist[i,1]=='X')
              {
                genos[genotypes$sex=='m' & genotypes$matRIL<21000,c("BB1","BB2","BB3","BB4","BB5","BB6","BB7","BB8")]<-NA
                genos[genotypes$sex=='m' & genotypes$matRIL>21000,c("AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8")]<-NA
              }
              
              big.list[[i]]<-genos
              rm(list=objnameA,pos=.GlobalEnv)
              rm(list=objnameB,pos=.GlobalEnv)
              
            }
            
          }else{
            stop("Design not valid. 
                 Must be one of inbredA, inbredB, ABcross, 
                 AAcross, or BBcross")
          } #last else
        }#AB if
      }#BB if
    }  #AA if
  }  #inbredB if
  
  
  
  
  #this makes sure all your phenotyped rils are in the matrix of genotypes
  phenotype.dat<-phenotype.dat[phenotype.dat$id %in% rownames(big.list[[1]]),]
  #order phenotype.dat by id
  phenotype.dat<-phenotype.dat[order(phenotype.dat$id),]
  return(list('genolist'=big.list,'positions','phenotype'=phenotype.dat))
    }#function close
