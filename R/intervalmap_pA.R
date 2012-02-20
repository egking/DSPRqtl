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

LocalInt<-(peakChr,peakPos,range=100,phenotype.dat,model,design)

ind.pos<-which(poslist, chr==peakChr & Ppos==peakPos)

pos.list<-poslist[(ind.pos-range):(ind.pos+range),]

output<-data.frame(
'chr'=numeric(nrow(pos.list)),
'Ppos'=numeric(nrow(pos.list)),
'Gpos'=numeric(nrow(pos.list)),
'LOD'=numeric(nrow(pos.list)))



null.mod<-lm(Median.Starv.ToD.h~1,data=phenotype)



for (i in 1:nrow(pos.list)) 

{	

filename<-paste("/home/eking/GenomeCache/A/",pos.list[i,1],"_",pos.list[i,2],".rda",sep="")
#filename<-paste("/home/eking/GenomeCache/A/2L_",14650000,".rda",sep="")

load(filename)


Agenotypes <-Ahmm_genos[,c("ril","AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8")]
patgeno <-merge(Agenotypes, phenotype,by.x="ril",by.y="patRIL")

X<-patgeno[,c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')]
Y<-patgeno[,'Median.Starv.ToD.h']

#get means for start
if(i ==1)
{
test<-summary(lm(patgeno$Median.Starv.ToD.h~1))
u1<-test$coef[1]
u2<-test$coef[1]
u3<-test$coef[1]
u4<-test$coef[1]
u5<-test$coef[1]
u6<-test$coef[1]
u7<-test$coef[1]
u8<-test$coef[1]

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

  resid0 <- lm(patgeno$Median.Starv.ToD.h~1)$resid  
  n.ind<-nrow(patgeno)
  nllik0 <- -sum(dnorm(resid0, 0, sqrt(sum(resid0^2)/n.ind), log=TRUE))/log(10)
cat(i,Sys.time(),nllik0-sLL.h1,"\n")

output[i,1]<-pos.list[i,1]
output[i,2]<-pos.list[i,2]
output[i,3]<-nllik0-sLL.h1

}

outname<-paste("/home/eking/QTL_results/AintervalSTARV_FEMALE",chr.all[j],"_",pos.all[j],".txt",sep="")
write.table(output,file=outname,row.names=FALSE)
cat(j,Sys.time(),"\n")
}



