##' \code{geno.means} estimates the founder genotype means and 
##' standard errors at a given position.
##' 
##' @title Founder genotype means and standard errors
##'   
##' @param peakChr character vector of length one. Must be one of the 
##'   major chromosome arms in the \emph{Drosophila} genome 
##'   ('X','2L','2R','3L',or '3R').
##'   
##' @param peakPos numeric vector of length one. A position in base 
##'   pairs in the DSPR position list (every 10kb).
##'   
##' @param model an object of class formula: a symbolic description of
##'   the null model to be fitted at each position (e.g., 
##'   \code{phenotype ~ 1}).  The genotype effects to be fitted will 
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
##' @param sex a character string (either 'm' or 'f') specifying the 
##'   sex of the measured individuals. This argument must be supplied 
##'   for a cross design for correct specification of the genotypes on
##'   the X chromosome.
##'   
##' @return A \code{data.frame} of the estimated mean and standard 
##'   error for each founder genotype. If a covariate is included in 
##'   the model statement, the estimate will be the founder genotype 
##'   mean after correcting for the covariate. If using the ABcross 
##'   design, a \code{list} of \code{data.frames} for pA and pB are 
##'   returned.
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
geno.means<-function(peakChr,peakPos,model,design,phenotype.dat, id.col,sex)
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
    
    patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
    mat<-as.matrix(patgeno[order(patgeno$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
    row.names(mat)<-patgeno$id
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
      
      patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
      mat<-as.matrix(patgeno[order(patgeno$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
      row.names(mat)<-patgeno$id
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
            mat<-as.matrix(matgeno[order(matgeno$id),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
            row.names(mat)<-matgeno$id
            rm(list=objnameA)
            rm(list=objnameB)
            
          }
          
          if(nrow(ABphenotype.dat)>0)
          {
            matgeno<-merge(ABphenotype.dat,Bgenotypes,by.x='matRIL',by.y='ril')
            matgeno<-merge(matgeno,Agenotypes,by.x='patRIL',by.y='ril') 
            mat<-as.matrix(matgeno[order(matgeno$id),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
            row.names(mat)<-matgeno$id
            rm(list=objnameA)
            rm(list=objnameB)
          }
          
          }else{           
            ABgenotypes<-merge(ABphenotype.dat,Agenotypes,by.x='patRIL',by.y='ril') 
            ABgenotypes<-merge(ABgenotypes, Bgenotypes, by.x='matRIL',by.y='ril',sort=FALSE)
            
            BAgenotypes<-merge(BAphenotype.dat,Agenotypes,by.x='matRIL',by.y='ril') 
            BAgenotypes<-merge(BAgenotypes, Bgenotypes, by.x='patRIL',by.y='ril',sort=FALSE)
            
            genotypes<-rbind(ABgenotypes,BAgenotypes)
            genotypes<-genotypes[order(genotypes$id),]
            
            mat<-as.matrix(genotypes[,c("AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8",
                                        "BB1","BB2","BB3","BB4","BB5","BB6","BB7","BB8")])
            row.names(mat)<-genotypes$id 
          }
      }else{
        stop("design must be one of: inbredA, inbredB, or ABcross")
      }#ab else  
        }#b else
    }#a else
  
  phenotype.dat<-phenotype.dat[order(phenotype.dat$id),]
  #this makes sure all your phenotyped rils are in the matrix of genotypes
  phenotype.dat<-phenotype.dat[phenotype.dat$id %in% rownames(mat),]
  
  #are there covariates?
  if(grepl("~\\s*1\\s*$", deparse(model)))
  {
  if(ncol(mat)==8){
    qtlmeanmod<-as.formula(paste(deparse(model),"+mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]+mat[,8]-1",sep=""))
    qtlmean<-lm(qtlmeanmod,data=phenotype.dat)
    qtlmean.sum<-summary(qtlmean)
    geno.means<-data.frame(qtlmean.sum$coef[(nrow(qtlmean.sum$coef)-7):nrow(qtlmean.sum$coef),1:2])   
    
    rownames(geno.means)<-colnames(mat)
    colnames(geno.means)<-c('Estimate','Std. Error')
    
    return(geno.means)
  }else{
    qtlmeanmod<-as.formula(paste(deparse(model),"+mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]+mat[,8]-1",sep=""))
    qtlmean<-lm(qtlmeanmod,data=phenotype.dat)
    qtlmean.sum<-summary(qtlmean)
    Ageno.means<-data.frame(qtlmean.sum$coef[(nrow(qtlmean.sum$coef)-7):nrow(qtlmean.sum$coef),1:2])   
    rownames(Ageno.means)<-colnames(mat)[1:8]
    colnames(Ageno.means)<-c('Estimate','Std. Error')
    
    qtlmeanmod<-as.formula(paste(deparse(model),"+mat[,9]+mat[,10]+mat[,11]+mat[,12]+mat[,13]+mat[,14]+mat[,15]+mat[,16]-1",sep=""))
    qtlmean<-lm(qtlmeanmod,data=phenotype.dat)
    qtlmean.sum<-summary(qtlmean)
    Bgeno.means<-data.frame(qtlmean.sum$coef[(nrow(qtlmean.sum$coef)-7):nrow(qtlmean.sum$coef),1:2])   
    rownames(Bgeno.means)<-colnames(mat)[9:16]
    colnames(Bgeno.means)<-c('Estimate','Std. Error')
    
    return(list(Ageno.means,Bgeno.means))
    
  }
  }else{
    #remove covariate effects
    p.resid<-lm(as.formula(model),data=phenotype.dat)$residuals
    pheno.name<-strsplit(deparse(model)," ~ ",fixed=TRUE)[[1]][1]    
    mean.pheno<-mean(phenotype.dat[,pheno.name])
    
    if(ncol(mat)==8){
      qtlmeanmod<-as.formula(p.resid~mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]+mat[,8]-1)
      qtlmean<-lm(qtlmeanmod,data=phenotype.dat)
      qtlmean.sum<-summary(qtlmean)
      geno.means<-data.frame(qtlmean.sum$coef[(nrow(qtlmean.sum$coef)-7):nrow(qtlmean.sum$coef),1:2])   
      
      rownames(geno.means)<-colnames(mat)
      colnames(geno.means)<-c('Estimate','Std. Error')
      geno.means['Estimate']<-geno.means[,'Estimate']+mean.pheno
      return(geno.means)
    }else{
      qtlmeanmod<-as.formula(p.resid~mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]+mat[,8]-1)
      qtlmean<-lm(qtlmeanmod,data=phenotype.dat)
      qtlmean.sum<-summary(qtlmean)
      Ageno.means<-data.frame(qtlmean.sum$coef[(nrow(qtlmean.sum$coef)-7):nrow(qtlmean.sum$coef),1:2])   
      rownames(Ageno.means)<-colnames(mat)[1:8]
      colnames(Ageno.means)<-c('Estimate','Std. Error')
      Ageno.means['Estimate']<-Ageno.means[,'Estimate']+mean.pheno
      
      qtlmeanmod<-as.formula(p.resid~mat[,9]+mat[,10]+mat[,11]+mat[,12]+mat[,13]+mat[,14]+mat[,15]+mat[,16]-1)
      qtlmean<-lm(qtlmeanmod,data=phenotype.dat)
      qtlmean.sum<-summary(qtlmean)
      Bgeno.means<-data.frame(qtlmean.sum$coef[(nrow(qtlmean.sum$coef)-7):nrow(qtlmean.sum$coef),1:2])   
      rownames(Bgeno.means)<-colnames(mat)[9:16]
      colnames(Bgeno.means)<-c('Estimate','Std. Error')
      Bgeno.means['Estimate']<-Bgeno.means[,'Estimate']+mean.pheno
      
      return(list(Ageno.means,Bgeno.means))
      
    }
  }#if/else covariate close
  }#function close