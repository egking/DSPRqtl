
#takes output from DSPRscan 

findQTL<-function(qtldat,threshold)
{
  output<-data.frame('chr'=numeric(length=0),'Ppos'=numeric(length=0),'Gpos'=numeric(length=0),'LOD'=numeric(length=0))
  for(i in 1:nrow(qtldat))
  {
    if(i==1 | i==nrow(qtldat)){}else{
      if(qtldat$LOD[i]>qtldat$LOD[i-1] & qtldat$LOD[i]>qtldat$LOD[i+1] & qtldat$LOD[i]>threshold)
      {
        peak<-data.frame('chr'=qtldat[i,'chr'],'Ppos'=qtldat[i,'Ppos'],'Gpos'=qtldat[i,'Gpos'],'LOD'=qtldat[i,'LOD'])
        output<-rbind(output,peak)
      }
    }
  }
  return(output)
}


findCI<-function(peakChr,peakPos,qtldat,peak,drop)
{
  startIndex<-which(qtldat$chr==peakChr & qtldat$Ppos==peakPos)
  l.ci<-peak-drop  #drop threshold
  start1<-startIndex #this one is for the other side of the peak
  #one side
  while(qtldat$LOD[start1]>l.ci & (start1-1) >= 1) 
  {  
    start1<-start1-1 
  }
  
  #other side
  while(qtldat$LOD[startIndex]>l.ci & (startIndex+1) <= nrow(qtldat))
  {  
    startIndex<-startIndex+1 
  }
  CIs<-qtldat[c(start1,startIndex),]
  rownames(CIs)<-c('Lower Bound','Upper Bound')
  return()
}

#means
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


founderNs<-function(peakChr,peakPos,model,design,phenotype.dat)
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





