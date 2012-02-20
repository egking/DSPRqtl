##' Effect size
##'
##' \code{perct.var} calculates the percent of the variance explained by the QTL.
##' @title Effect size
##' @param peakChr character vector of length one. 
##' Must be one of the major chromosome arms in the 
##' \emph{Drosophila} genome ('X','2L','2R','3L',or '3R').
##' 
##' @param peakPos numeric vector of length one. A position in 
##' base pairs in the DSPR position list (every 10kb). See 
##' \code{\link{positionlist_wgenetic}}.
##' 
##' @param model an object of class formula: a symbolic description of the 
##' null model to be fitted at each position (e.g., \code{phenotype ~ 1}).  The genotype effects to be fitted 
##' will be added based on \code{design}. 
##' 
##' @param phenotype.dat \code{data.frame} containing a column of 
##' ril ids (must be named patRIL) and phenotypes.
##' 
##' @param design a character string. One of either 'inbredA' or 'inbredB'
##' corresponding to the pA and pB set of inbred RILs. Other crossing designs 
##' will be supported in the future.
##' 
##' @return A numeric vector of length one: the percent variance explained by the QTL.
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##' 
perct.var<-function(peakChr,peakPos,model,phenotype.dat,design)
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
      
      patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
      mat<-as.matrix(patgeno[order(patgeno$patRIL),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
      row.names(mat)<-patgeno$patRIL
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
      
      patgeno<-merge(phenotype.dat,genotypes,by.x="patRIL",by.y="ril")
      mat<-as.matrix(patgeno[order(patgeno$patRIL),c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')])
      row.names(mat)<-patgeno$patRIL
    }
    
  }
  phenotype.dat<-phenotype.dat[order(phenotype.dat$patRIL),]
  #this makes sure all your phenotyped rils are in the matrix of genotypes
  phenotype.dat<-phenotype.dat[phenotype.dat$patRIL %in% rownames(mat),]
  
  qtlmodel<-as.formula(paste(deparse(model),"+mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]",sep=""))
  qtlmod<-lm(qtlmodel,data=phenotype.dat)
  sumqtl<-summary.aov(qtlmod)
  sumqtl<-as.data.frame(sumqtl[[1]])
  ss<-sumqtl["mat[, 1]","Sum Sq"]+sumqtl["mat[, 2]","Sum Sq"]+sumqtl["mat[, 3]","Sum Sq"]+sumqtl["mat[, 4]","Sum Sq"]+
    sumqtl["mat[, 5]","Sum Sq"]+sumqtl["mat[, 6]","Sum Sq"]+sumqtl["mat[, 7]","Sum Sq"]
  
  ndf<-7
  ddf<-sumqtl["Residuals","Df"]
  mss<-ss/ndf
  mse<-sumqtl["Residuals","Mean Sq"]
  test.f <-mss/mse
  return(100*(test.f/(test.f+(ddf/ndf))))
}
