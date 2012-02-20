load(file='/home/eking/QTL_results/GD_prelim_rank.rda')
#q1<-lod.matrix[[3]]
#load(file='/home/eking/QTL_results/GD_prelim_raw.rda')
#q2<-lod.matrix[[3]]


positions<-read.table(file="/home/eking/QTL_results/Starv_mA_20Dec_RESULTS.txt",header=TRUE)
positions<-positions[,c('Gaxis')]

#READ IN DATA (CHANGE PATH TO YOUR FILE)

qtl.data<-lod.matrix[[1]]



##MAKE GENOME-WIDE PLOT

#size of the plot - use pdf() to write out to a file:: might look a little better if it is a bit wider
#pdf("/Users/elizabethking/Dropbox/SJM_EGK/QTL_Data/Nicotine/NicotineA.pdf",height = 5, width=12)
pdf("/home/eking/QTL_results/D1_rank.pdf",height = 5, width=12)

#quartz(width=8,height=5)
#mar changes the size of the margins, mpg changes the placement of the axis labels
par(mar=c(7,7,4,2),mgp=c(4.5,1,0))

#when plotting more than one (e.g. males and females on same plot), 
#use the one with the higher LOD to make the y scale correctly
plot(positions,qtl.data,
	xlab="Position (cM)",
	ylab="LOD",
	cex.lab=1.5, #size of text
	type="n",	#type = none: don't make plot yet, need to make rectangles first
	axes=FALSE  #we will do the axes by hand next
)

#make grey boxes to delineate the different chromosome arms
	#grey90 is only 10% grey- decrease the number to make it darker
	#must be placed before lines or they will cover the lines
	#again, use the data with the highest LOD
top<-max(qtl.data)

rect(66.3,0,120.3,top,col="grey90",lty=0)
rect(174,0,221,top,col="grey90",lty=0)

points(positions,qtl.data,
		type="l",
		lwd=2, #thickness of line
		col="black")		
#* points(qtl.data2$Gaxis,qtl.data2$LOD,
#*		type="l",
#*		lwd=2,
#*		col="red")		
			
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
#add horizontal line at LOD of 4
abline(h=10.1)

#what colors did you use?
#* colors<-c("black","red")
#* legend("topleft",legend=c("females","males"),col=colors,lwd=2,cex=1.5)
dev.off()