##' \code{LocalInt} performs standard interval mapping around a given 
##' peak for a range of positions specified by the user.
##' 
##' @title Local Interval Mapping
##'   
##' @param peakChr character vector of length one. Must be one of the 
##'   major chromosome arms in the \emph{Drosophila} genome 
##'   ('X','2L','2R','3L',or '3R').
##'   
##' @param peakPos numeric vector of length one. A position in base 
##'   pairs in the DSPR position list (every 10kb).
##'   
##' @param range numeric vector of length one specifying the number of
##'   positions to test on either side of the peak (positions are 
##'   every 10kb). Default is 100.
##'   
##' @param phenotype.dat \code{data.frame} containing a column of ril 
##'   ids (must be named patRIL) and phenotypes.
##'   
##' @param pheno.name a character string specifying the column name of
##'   the phenotype data
##'   
##' @param design a character string. One of either 'inbredA' or 
##'   'inbredB' corresponding to the pA and pB set of inbred RILs. 
##'   Other crossing designs may be supported in the future.
##'   
##' @return A \code{data.frame} consisting of the chromosome, physical
##'   position (bp), genetic position (cM) and LOD score for each 
##'   position.
##' 
##' @author Elizabeth King (\email{kingeg@@missouri.edu})
##' 
##' 
##' @export
##' 
LocalInt<-function(peakChr,peakPos,range=100,phenotype.dat,pheno.name,design)
{
  
  LogLik<-function(x,p)
  {
    uis<-p[1:8]
    sigma<-p[9]
    Xis<-x[1:8]
    Yi<-unlist(rep(x[9],8))
    d.stat<-Yi-uis
    l.vect<-Xis*dnorm(d.stat,sd=sigma)
    if(sum(l.vect==0)){l.vect<-0.00001}
    return(log10(sum(l.vect)))
    
  }
  
  SumLogLik<-function(p)
  {
    sLL<- -sum(apply(cbind(X,Y),1,function(x) LogLik(x,p)))
    return(sLL)
  }
  
  
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

  data(positionlist_wgenetic,envir=environment())
  names(poslist)<-c('chr','Ppos','Gpos','Gaxis')
  
  ind.pos<-which(poslist$chr==peakChr & poslist$Ppos==peakPos)
  
  int.list<-poslist[(ind.pos-range):(ind.pos+range),]
  
  output<-data.frame(
    'chr'=numeric(nrow(int.list)),
    'Ppos'=numeric(nrow(int.list)),
    'Gpos'=numeric(nrow(int.list)),
    'LOD'=numeric(nrow(int.list)))

 
  for (i in 1:nrow(int.list)) 
    
  {  
    

  if(design=='inbredA'|design=='inbredB')
  {
    
    if(design=='inbredA')
    {
      
      objname<-paste("A_",int.list[i,1],"_",format(int.list[i,2], sci = FALSE),sep="")
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
      X<-patgeno[,c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')]
      Y<-patgeno[,pheno.name]
      rm(list=objname)
     
    }else{
      objname<-paste("B_",int.list[i,1],"_",format(int.list[i,2], sci = FALSE),sep="")
      
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
      X<-patgeno[,c('BB1','BB2','BB3','BB4','BB5','BB6','BB7','BB8')]
      Y<-patgeno[,pheno.name]
      rm(list=objname)
      
    }
    
  }
  
nmod<-Y~1

#get means for start
if(i ==1)
{
test<-summary(lm(Y~X[,1]+X[,2]+X[,3]+X[,4]+X[,5]+X[,6]+X[,7]+X[,8]-1))
u1<-test$coef[1]
u2<-test$coef[2]
u3<-test$coef[3]
u4<-test$coef[4]
u5<-test$coef[5]
u6<-test$coef[6]
u7<-test$coef[7]
u8<-test$coef[8]

sigma<-test$sigma
}else{
 u1<-p.new[1]
  u2<-p.new[2]
  u3<-p.new[3]
  u4<-p.new[4]
  u5<-p.new[5]
  u6<-p.new[6]
  u7<-p.new[7]
  u8<-p.new[8]
  sigma<-p.new[9]
}
  
#optim
out<-optim(c(u1,u2,u3,u4,u5,u6,u7,u8,sigma),SumLogLik,method="L-BFGS-B",lower=c(rep(mean(Y)-2*sd(Y),8),1e-2),upper=c(rep(mean(Y)+2*sd(Y),8),sd(Y)))

p.new<-c(out$par[1],out$par[2],out$par[3],out$par[4],out$par[5],out$par[6],out$par[7],out$par[8],out$par[9])

sLL.h0<-SumLogLik(c(rep(mean(Y),8),p.new[9]))
sLL.h1<-SumLogLik(p.new)

  
  resid0 <- lm(nmod)$resid  
  n.ind<-nrow(patgeno)
  nllik0 <- -sum(dnorm(resid0, 0, sqrt(sum(resid0^2)/n.ind), log=TRUE))/log(10)

output[i,1]<-int.list[i,1]
output[i,2]<-int.list[i,2]
output[i,3]<-int.list[i,3]
output[i,4]<-nllik0-sLL.h1

}#for close
return(output)
}#function close


