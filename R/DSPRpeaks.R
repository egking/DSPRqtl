DSPRpeaks<-function(qtldat,threshold,LODdrop)
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
    Ns<-founderNs(peakChr,peakPos,model,design)
    entropy<-entropy.pos(peakChr,peakPos,phenotype.dat)
    peaks.obj[[p]]<-list(
                    "threshold"=threshold,
                    "peak"=peakmat[p,],
                    "LODdrop"=LODdrop
                    "CI"=ci,
                    "founderNs"=Ns,
                    "geno.means"=gmeans,
                    "perct.var"=pvar,
                    "entropy"=entropy                
                    )
  }
  
class(peaks.obj)<-'peaks'
return(peaks.obj)
  
}



print.peaks <- function(x, digits = 4, ...){
  cat("\n")
  cat("Significant peaks at a threshold of",x[[1]]$threshold,"\n")
  cat("\tAll peaks should be examined visually by the user to determine which constitute distinct peaks.\n")
  cat("\n")
  for (i in 1:length(x))
  {
  xi<-x[[i]]  
  
  cat("Peak ",i,":\n")
  cat("\t",xi$peak,"\n")
  cat("\n")
  
  cat("Confidence Interval(LOD drop = ",xi$LODdrop,"):\n")
  cat("\t",xi$CI,"\n")
  cat("\n")
  
  cat("Numbers of RILs with each founder genotype:\n")
  cat("\t",xi$Ns,"\n")
  cat("\n")
  
  cat("Founder genotype means and standard errors:\n")
  cat("\t",xi$geno.means,"\n")
  cat("\n")
  
  cat("Percent of variation explained by the QTL:\n")
  cat("\t",xi$perct.var,"\n")
  cat("\n")
  
  cat("Proportion of missing information (entropy):\n")
  cat("\t",xi$entropy,"\n")
  cat("\n\n\n")
}
}