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
##'   signficance threshold. Default is 6.8. Use
##'   \code{\link{DSPRperm}} to get a threshold specific to a given
##'   dataset.
##'   
##' @param LODdrop numeric vector of length one consisting of the LOD
##'   drop to be used. Default value is 2 which approximates a 95\%
##'   confidence interval.
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
DSPRpeaks<-function(qtldat,threshold=6.8,LODdrop=2)
{
  peakmat<-findQTL(qtldat$LODscores,threshold)
  peaks.obj<-vector('list',nrow(peakmat))
  for (p in 1:nrow(peakmat))
  {
    peak<-peakmat[p,'LOD']
    peakChr<-peakmat[p,'chr']
    peakPos<-peakmat[p,'Ppos']
    model<-qtldat$model
    design<-qtldat$design
    phenotype.dat<-qtldat$phenotype
    ci<-findCI(peakChr,peakPos,qtldat$LODscores,peak,LODdrop)
    gmeans<-geno.means(peakChr,peakPos,model,phenotype.dat,design)
    pvar<-perct.var(peakChr,peakPos,model,phenotype.dat,design)
    Ns<-founderNs(peakChr,peakPos,design,phenotype.dat)
    ent<-entropy.pos(peakChr,peakPos,phenotype.dat,design)
    peaks.obj[[p]]<-list(
                    "threshold"=threshold,
                    "peak"=peakmat[p,],
                    "LODdrop"=LODdrop,
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
  
  cat("Confidence Interval(LOD drop = ",xi$LODdrop,"):\n")
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
