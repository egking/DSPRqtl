##' \code{DSPRpeaks} takes output from \code{\link{DSPRscan}}. Locates
##' and summarizes QTL peaks.
##' 
##' @title Find and Summarize QTL
##'   
##' @aliases DSPRpeaks print.peaks
##'   
##' @param qtldat An object of class gscan. Output from
##'   \code{\link{DSPRscan}}.
##'   
##' @param threshold numeric vector of length one consisting of the 
##'   signficance threshold. Default is 6.8 for inbred designs and 10.1 
##'   for the ABcross. Use \code{\link{DSPRperm}} to get a threshold 
##'   specific to a given dataset.
##'
##' @param method a character string specifying the method to use for 
##'   the confidence interval. Options are: 'LODdrop' calculates a LOD 
##'   drop interval for the drop amount specified, 'BCI' calculates the 
##'   Bayesian credible interval for the fraction specified, and 'both' 
##'   calculates and returns both. 
##'      
##' @param LODdrop numeric vector of length one consisting of the LOD
##'   drop to be used when using method 'LODdrop'. Default value is 2 
##'   which approximates a 95\% confidence interval for the inbred 
##'   designs. Users of the ABcross design should consider using a 
##'   larger LOD drop.
##'   
##' @param BCIprob numeric vector of length one consisting of the nominal
##' Bayes fraction. Default value is 0.95. 
##' 
##' @return A list of class \code{peaks} containing a list for each
##'   significant peak, each containing:
##' \item{threshold}{the specified signficance threshold}
##' \item{peak}{A single row \code{data.frame} with the chromosome, 
##' physical position (bp), genetic position (cM) and LOD score for
##' the peak}
##' \item{LODdrop}{the specified LOD drop for the support interval}
##' \item{CI}{the upper and lower bounds of the confidence interval}
##' \item{founderNs}{the number of RILs with each founder genotype at the peak}
##' \item{geno.means}{the estimated means and standard errors for each founder genotype}
##' \item{perct.var}{the percent variation explained by the QTL}
##' \item{entropy}{the proporiton of missing information (entropy)}
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##'
##' @S3method print peaks
DSPRpeaks<-function(qtldat,method,threshold,LODdrop,BCIprob)
{
  if(missing(method))
  {
    method<-'both'
  }
  
  if(missing(LODdrop))
  {
    if(method=='LODdrop'|method=='both')
    {
      LODdrop<-2
    }else{
      LODdrop<-NA
    } 
  }  
  
  if(missing(BCIprob))
  {
    if(method=='BCI'|method=='both')
    {
      BCIprob<-0.95
    }else{
      BCIprob<-NA
    } 
  }  
  
  if(missing(threshold))
  {
    if(qtldat$design=='ABcross')
    {
      threshold<-10.1
    }else{
      threshold<-6.8
    }
  }
  
  peakmat<-findQTL(qtldat$LODscores,threshold)
  peaks.obj<-vector('list',nrow(peakmat))
  for (p in 1:nrow(peakmat))
  {
    peakChr<-peakmat[p,'chr']
    peakPos<-peakmat[p,'Ppos']
    model<-qtldat$model
    design<-qtldat$design
    phenotype.dat<-qtldat$phenotype
    ci<-findCI(peakChr,peakPos,qtldat$LODscores,method,LODdrop,BCIprob)
    gmeans<-geno.means(peakChr,peakPos,model,design,phenotype.dat,id.col='id',sex=qtldat$sex)
    pvar<-perct.var(peakChr,peakPos,model,design,phenotype.dat,id.col='id',sex=qtldat$sex)
    Ns<-founderNs(peakChr,peakPos,design,phenotype.dat,id.col='id')
    ent<-entropy.pos(peakChr,peakPos,design,phenotype.dat,id.col='id',sex=qtldat$sex)
    peaks.obj[[p]]<-list(
                    "threshold"=threshold,
                    "peak"=peakmat[p,],
                    "LODdrop"=LODdrop,
                    "BCIprob"=BCIprob,
                    "CI"=ci,
                    "founderNs"=Ns,
                    "geno.means"=gmeans,
                    "perct.var"=pvar,
                    "entropy"=ent                
                    )
  }
  
class(peaks.obj)<-'peaks'
return(peaks.obj)
  
}



print.peaks <- function(x, ...){
  cat("\n")
  cat("Significant peaks at a threshold of",x[[1]]$threshold,"\n")
  cat("\tAll peaks should be examined visually by the user to determine which constitute distinct peaks.\n")
  cat("\n")
  for (i in 1:length(x))
  {
  xi<-x[[i]]  
  
  cat("Peak ",i,":\n")
  print(xi$peak)
  cat("\n")
  
  cat("Confidence Interval(LOD drop = ",xi$LODdrop,", BCI fraction = ",xi$BCIprob,"):\n")
  print(xi$CI)
  cat("\n")
  
  cat("Numbers of RILs with each founder genotype:\n")
  print(xi$founderNs)
  cat("\n")
  
  cat("Founder genotype means and standard errors:\n")
  print(xi$geno.means)
  cat("\n")
  
  cat("Percent of variation explained by the QTL:\n")
  cat("\t",xi$perct.var,"\n")
  cat("\n")
  
  cat("Proportion of missing information (entropy):\n")
  cat("\t",xi$entropy,"\n")
  cat("\n\n\n")
}
}
