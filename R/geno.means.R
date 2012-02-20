geno.means<-function(peakChr,peakPos,model,phenotype.dat,design)
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
  
  qtlmeanmod<-as.formula(paste(deparse(model),"+mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]+mat[,8]-1",sep=""))
  qtlmean<-lm(qtlmeanmod,data=phenotype.dat)
  qtlmean.sum<-summary(qtlmean)
  geno.means<-data.frame(qtlmean.sum$coef[,1:2])
  if(design=='inbredA'|design=='inbredB')
  {
    
    if(design=='inbredA')
    {
      rownames(geno.means)<-c('A1','A2','A3','A4','A5','A6','A7','A8')
      colnames(geno.means)<-c('Estimate','Std. Error')
    }else{
      rownames(geno.means)<-c('B1','B2','B3','B4','B5','B6','B7','B8')
      colnames(geno.means)<-c('Estimate','Std. Error')
    }
    
  }
  return(geno.means)
}
