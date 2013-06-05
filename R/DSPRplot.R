##' Function to plot genome scan results for the DSPR RILs
##'
##' @title Genome Scan Plot
##' 
##' @param qtldata a \code{list} of output from DSPRscan with each 
##'   list element corrsponding to one DSPRscan result. e.g. to plot 
##'   multple genome scans on a single plot: \code{qtldata <-
##'   list(scan1,scan2)}
##'   
##' @param threshold numeric vector of length one consisting of the 
##'   signficance threshold. Default is 6.8 for inbred designs and 10.1 
##'   for the ABcross. Use \code{\link{DSPRperm}} to get a threshold 
##'   specific to a given dataset.
##'   dataset.
##'   
##' @param legNames a character vector with names for each DSPRscan 
##'   result to be plotted. Defaults to the phenotype column names
##'   used in DSPRscan.
##'   
##' 
##' @author Elizabeth King (\email{egking@@uci.edu})
##' 
##' @export
##'
DSPRplot<-function(qtldata,threshold,legNames=NULL)
{
  
  if(missing(threshold))
  {
    if(qtldat$design=='ABcross')
    {
      threshold<-10.1
    }else{
      threshold<-6.8
    }
  }
  
  
  
if(length(names(qtldata))>0){qtldata<-list(qtldata)}
  
dots <- list(...)

data(positionlist_wgenetic)

LODscores<-vector("list",length=length(qtldata))
allLODs<-numeric(0)
if(is.null(legNames)){
  legNames<-numeric(0)
  legtest<-TRUE
  }
for(i in 1:length(qtldata))
{
LODscores[[i]]<-merge(qtldata[[i]]$LODscores,poslist, by=c('chr','Ppos','Gpos'),sort=FALSE)
allLODs<-c(allLODs,qtldata[[i]]$LODscores$LOD)
if(legtest){legNames<-c(legNames,as.character(qtldata[[i]]$model)[2])}
}

top<-max(allLODs)

allcols<-rep(c('black','red','blue','green','orange','purple','cyan','brown'),5)

#quartz(width=5,height=12)
par(mar=c(7,7,4,2),mgp=c(4.5,1,0))

plot(LODscores[[1]]$Gaxis, LODscores[[1]]$LOD,
	xlab="Position (cM)",
	ylab="LOD",
  ylim=c(0,top),
	cex.lab=1.5, #size of text
	type="n",	#type = none: don't make plot yet, need to make rectangles first
	axes=FALSE  #we will do the axes by hand next
)

rect(66.3,-10,120.3,top+10,col="grey90",lty=0)
rect(174,-10,221,top+10,col="grey90",lty=0)

for(j in 1:length(LODscores))
{
points(LODscores[[j]]$Gaxis, LODscores[[j]]$LOD,
		type="l",
		lwd=2, #thickness of line
		col=allcols[j])		
}			
box() #add box back in

axis(1, at = c(0,66.3,120,174,221,277), labels= c(0,"66  0",54,"108  0",47,103), tick=FALSE,cex.axis=1)

# add y axis back in, cex.axis=size of numbers 
axis(2,cex.axis=1.5)

#add chromosome labels at midpts
mtext("X",line = 2.5,side=1, at =33.15, cex=1.5)
mtext("2L",line = 2.5,side=1, at =93.3, cex=1.5)
mtext("2R",line = 2.5,side=1, at =147, cex=1.5)
mtext("3L",line = 2.5,side=1, at =197.5, cex=1.5)
mtext("3R",line = 2.5,side=1, at =249, cex=1.5)
#add horizontal line at LOD of 6.8
abline(h=threshold)

#what colors did you use?
legcol<-allcols[1:length(LODscores)]
if(length(LODscores)>1){legend("topleft",legend=legNames,col=legcol,lwd=2,cex=1.1)}
}
