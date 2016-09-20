setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model")

load("Farm model/farm_fits_20160209.RData")

library(gplots)
library(dclone)


T.all<-as.Date(c(min(lice_data$date):max(lice_data$date)), origin="1970-01-01") # Note that origin may be different on different operating systems (set for Mac)

p.fitted<-list(); length(p.fitted)<-3
sum.fits<-list(); length(sum.fits)<-3
for(i in 1:3){
	sum.fits[[i]]<-summary(allFits[[i]])[[1]]
	p.fitted[[i]]<-matrix(NA, nrow=4, ncol=2)
	rownames(p.fitted[[i]])<-rownames(sum.fits[[i]])
	colnames(p.fitted[[i]])<-c("mean", "sd")
	for(j in 1:4){
		p.fitted[[i]][j,1]<-sum.fits[[i]][j,1]
		p.fitted[[i]][j,2]<-sum.fits[[i]][j,3]
		}
}

y.pred<-list(); length(y.pred)<-3
dt<-1 #timestep in simulations is 1 day
for(i in 1:3){
	y.pred[[i]]<-matrix(nrow=length(T.all), ncol=3)
	ind1<-which(T.all<treatment.date[i])
	ind2<-which(T.all>=treatment.date[i])
	
	r1<-rep(p.fitted[[i]][3,1], length(ind1))*dt
	r1.var<-rep((p.fitted[[i]][3,2])^2, length(ind1))*dt
	r2<-rep(p.fitted[[i]][4,1], length(ind2))*dt
	r2.var<-rep((p.fitted[[i]][4,2])^2, length(ind2))*dt
	
	y.pred[[i]][ind1,1]<-p.fitted[[i]][1,1]*exp(r1*as.numeric(T.all[ind1]-treatment.date[i]))
	y.pred[[i]][ind2,1]<-p.fitted[[i]][1,1]*exp(r2*as.numeric(T.all[ind2]-treatment.date[i]))
	
	#lower:
	y.pred[[i]][ind1,2]<-(p.fitted[[i]][1,1]-1.96*p.fitted[[i]][1,2])*exp((r1+1.96*sqrt(r1.var))*as.numeric(T.all[ind1]-treatment.date[i]))
	y.pred[[i]][ind2,2]<-(p.fitted[[i]][1,1]-1.96*p.fitted[[i]][1,2])*exp((r2-1.96*sqrt(r2.var))*as.numeric(T.all[ind2]-treatment.date[i]))
	
	#upper:
	y.pred[[i]][ind1,3]<-(p.fitted[[i]][1,1]+1.96*p.fitted[[i]][1,2])*exp((r1-1.96*sqrt(r1.var))*as.numeric(T.all[ind1]-treatment.date[i]))
	y.pred[[i]][ind2,3]<-(p.fitted[[i]][1,1]+1.96*p.fitted[[i]][1,2])*exp((r2+1.96*sqrt(r2.var))*as.numeric(T.all[ind2]-treatment.date[i]))
	
}

boot<-function(x, n=1000){
	xi<-matrix(sample(x, n*length(x), replace=TRUE), nrow=n, ncol=length(x))
	mean.x<-apply(xi, 1, mean, na.rm=TRUE)
	ret<-c(mean(x, na.rm=TRUE), quantile(mean.x, c(0.025, 0.975), na.rm=TRUE))
	return(ret)
}

y.true<-list(); length(y.true)<-3
for(i in 1:2){
	farm.ind<-which(lice_data$farm==c("SP", "HR", "BG")[i])
	sampdate<-unique(as.numeric(lice_data$date[farm.ind]))
	
	y.true[[i]]<-matrix(nrow=length(sampdate), ncol=4)
	y.true[[i]][,1]<-sampdate
	y.true[[i]][,2:4]<-matrix(unlist(tapply(lice_data$lice[farm.ind], as.numeric(lice_data$date[farm.ind]), boot)), length(sampdate), 3, byrow=TRUE) 		
}
farm.ind<-which(lice_data$farm=="BG")
y.true[[3]]<-cbind(lice_data$date[farm.ind], lice_data$lice[farm.ind]/20)
		
quartz(width=4, height=7)
par(mfrow=c(3,1), mar=c(4,3,1,1), oma=c(2,2,1,1))

for(i in 1:3){
	plot(T.all, y.pred[[i]][,1], "l", xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2)
	polygon(x=c(T.all, rev(T.all)), y=c(y.pred[[i]][,2], rev(y.pred[[i]][,3])), border=NA, col="#00000030")
	if(i<3) plotCI(y.true[[i]][,1], y.true[[i]][,2], li=y.true[[i]][,3], ui=y.true[[i]][,4], add=TRUE, gap=0.3, pch=21, pt.bg="white")
	if(i==3) points(y.true[[i]][,1], y.true[[i]][,2], pch=21, bg="white")
	abline(v=treatment.date[i], lty=2)
	}

mtext(side=1, outer=TRUE, "Date (2005/2006)")
mtext(side=2, outer=TRUE, expression(paste("Average motile ",  italic(L.), " ", italic(salmonis),  " per farmed salmon", sep="")))
	
	
#######################################################
# Data cloning results
#######################################################
DC<-list();length(DC)<-3
for(i in 1:3) DC[[i]]<-dctable(allFits[[i]])

par.names<-c(expression(f[0]), "k", expression(r[1]), expression(r[2]))
farm.name<-c("Sargeaunts", "Humphrey", "Burdwood")
quartz(width=6.3, height=5, pointsize=11)
par(mfrow=c(3,4), mar=c(3,3,0,0), oma=c(2,4,3,1))
for(j in 1:3){
	for(i in 1:4){
		plot(DC[[j]][[i]]$n.clones, DC[[j]][[i]]$sd^2/(DC[[j]][[i]]$sd[1]^2), "l", bty="l", las=1, ylab="", xlab="", las=1)
		lines(1:10, 1/(1:10), lty=2)
		points(DC[[j]][[i]]$n.clones, DC[[j]][[i]]$sd^2/(DC[[j]][[i]]$sd[1]^2), pch=21, bg="white")
		if(j==1) mtext(side=3, par.names[i], line=1)
		if(i==1) mtext(side=2, line=5, farm.name[j])
		if(i==1&j==2) mtext(side=2, line=3, "Scaled variance")
	}
}
mtext(side=1, outer=TRUE, "Number of clones", line=0)
segments(-42.5,0,-42.5,1.3,xpd=NA)
segments(-42.5,2.6,-42.5,3.9,xpd=NA)