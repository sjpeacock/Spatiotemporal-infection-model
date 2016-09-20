rm(list=ls())
setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model")
library(gplots)
#################################################################################################
## Set up parameters for treatment scenarios
#################################################################################################

source('params.R')
p<-p.lice[1,]


source("Simulations/sim-grid.R")
source("Simulations/sim-functions.R")
source("Fitting/sim-model.R")

# time zero is September 1, 2005 for our purposes
r1<-c(p.farm[[1]][1,1], p.farm[[2]][1,1], p.farm[[3]][1,1])
r2<-c(p.farm[[1]][2,1], p.farm[[2]][2,1], p.farm[[3]][2,1])

#-----------------------------------------------
# Scenarios are:
fT<-list(); length(fT)<-4
tT<-list(); length(tT)<-4

# 1) true treatment dates and abundances
fT[[1]]<-c(p.farm[[1]][3,1], p.farm[[2]][3,1], p.farm[[3]][3,1])
tT[[1]]<-c(as.numeric(as.Date("2006-03-07")), as.numeric(as.Date("2006-02-13")), as.numeric(as.Date("2006-01-13")))-as.numeric(as.Date("2005-09-01"))

# 2) treated independently when threshold crossed

#What was the predicted abundance on Sept 1?  This is where we start simulations
f0<-fT[[1]]*exp(-r1*tT[[1]])
fT[[2]]<-c(3,3,3)
tT[[2]]<-1/r1*log(fT[[2]]/f0)

# 3) all treated when the first farm crossed threhsold (information sharing)  
tT[[3]]<-rep(min(tT[[2]]), 3)
fT[[3]]<-f0*exp(r1*tT[[3]])

# 4) all treated Feb. 1
tT[[4]]<-rep((as.numeric(as.Date("2006-02-01"))-as.numeric(as.Date("2005-09-01"))), 3)
fT[[4]]<-f0*exp(r1*tT[[4]])

#-----------------------------------------------
# Plot louse abundances under the four scenarios

times<-c(0:303)
predicted.f<-list(); length(predicted.f)<-4
for(i in 1:4){
	predicted.f[[i]]<-matrix(NA, nrow=length(times), ncol=3)
	for(k in 1:3){
		for(j in 1:length(times)){
		predicted.f[[i]][j,k]<-fset(tT[[i]][k], fT[[i]][k], r1[k], r2[k], times[j])
		}
	}
}

quartz(width=5.3, height=4, pointsize=10)
par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(0,1,0,0))
topText<-c("a) Scenario A", "b) Scenario B", "c) Scenario C", "d) Scenario D")
for(i in 1:4){
	plot(as.Date(times, origin="2005-09-01"), predicted.f[[i]][,1], "l",col=1, ylim=c(0, 8), bty="l", las=1, xlab="", ylab="")#max(predicted.f[[i]])
	polygon(x=c(as.Date("2006-04-01"), as.Date("2006-04-01"), as.Date("2006-06-30"), as.Date("2006-06-30")), y=c(-1, 10, 10, -1), border=NA, col="#00000020")
	lines(as.Date(times, origin="2005-09-01"), predicted.f[[i]][,2], lty=2)
	lines(as.Date(times, origin="2005-09-01"), predicted.f[[i]][,3], lty=3)
	#abline(v=as.Date(tT[[i]], origin="2005-09-01"), lty=c(1:3))
	abline(h=3, col=3)#, col=grey(0.8), lwd=2)
	mtext(side=3, adj=0, line=0.5, topText[i])
	if(i==2) legend("topleft", lty=c(1,2,3), c("Sargeaunts", "Humphrey", "Burdwood"), bg="white", bty="n")#legend(as.Date("2005-09-01"), 13, lwd=c(1,1,1,1,15), col=c(rep(1,3), 2,grey(0.8)), lty=c(1,2,3,1,1), c("Sargeaunts", "Humphrey", "Burdwood", "Treatment threshold", "Juvenile salmon migration"), bg="white", ncol=3, xpd=NA, bty="n")
}
mtext(side=2, outer=TRUE, expression(paste("Average number of motile ", italic(L), ". ",  italic(salmonis), " per farmed salmon")), line=-0.5)
#################################################################################################
## Create cope distributions for all scenarios
#################################################################################################
source('~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model/Simulations/sim-larvae.R')

migDate<-function(x, Xstart=-20, Tstart=195, v=p.lice[1,5], numeric=FALSE){
	y<-Tstart+(x-Xstart)/v
	if(numeric==FALSE) return(as.Date(y, origin="2006-01-01")) else return(y)
}

#Average dates of emergence (2007-2012, and 2010)
e1<-98
e2<-87

# copedistributions saved as list of length 4 (4 scenarios) called farmL.all
source('~/Desktop/filledContour.R')
pal<-colorRampPalette(c("white", 1, 1))
quartz(width=6.3, height=6, pointsize=12)
layout(matrix(c(1,3,1,3,1,3,2,4,2,4,2,4,5,5), 2, 7))
par(mar=c(3,3,2,1), oma=c(2,2,1,0))
for(i in 1:4){
	filledContour(farmL.all[[i]], x=x, y=as.Date(T, origin="2005-09-01"), zlim=c(0,0.05), color.palette=pal, xlab="", ylab="")
	contour(farmL.all[[i]], x=x, y=as.Date(T, origin="2005-09-01"), add=TRUE)
	mtext(side=3, line=0.05, paste("Scenario ", LETTERS[i]))
	arrows(x0=-40, y0=as.Date(paste("2006-", e1, sep=""), format="%Y-%j"), x1=60, y1=migDate(60, Xstart=-40, Tstart=e1), col=4, length=0.08, xpd=NA)
		points(-40, as.Date(paste("2006-", e1, sep=""), format="%Y-%j"), pch=21, col=4, bg="white")
	arrows(x0=-40, y0=as.Date(paste("2006-", e2, sep=""), format="%Y-%j"), x1=60, y1=migDate(60, Xstart=-40, Tstart=e2), col=2, length=0.08, xpd=NA)
	points(-40, as.Date(paste("2006-", e2, sep=""), format="%Y-%j"), pch=21, col=2, bg="white")
	
	
	if(i>=3) abline(h=as.Date(tT[[i]][1], origin="2005-09-01"), lty=2)
	abline(v=c(0,3.5,53), lty=c(1:3), lwd=1.2)	
	
	}
mtext(side=1, outer=TRUE, "Distance along migration (km)")
mtext(side=2, outer=TRUE, "Date")

plot.scale.bar(zlim=c(0, 0.05), color.palette=pal)
mtext(side=4, "Larval density", line=3.5)

#################################################################################################
## Calculate metrics
#################################################################################################
source('Fitting/sim-model.R')

set.seed(98375)

G<-read.csv("Simulations/Glendale.csv")
D1<-sample(G$day[G$year!=2010], size=1000, replace=TRUE, prob=G$fry[G$year!=2010]/sum(G$fry[G$year!=2010])) 
D2<-sample(G$day[G$year==2010], size=1000, replace=TRUE, prob=G$fry[G$year==2010]/sum(G$fry[G$year==2010])) 

quartz(width=3.2, height=2.5, pointsize=10)
hist(D1, col="#0000FF80", border="white", main="", xaxt="n", las=1, xlab="Date of emergence")
axis(side=1, at=seq(70, 120, 10), labels= strftime(as.Date(seq(70, 120, 10), origin="2006-01-01"), format="%b-%d"))
hist(D2, col="#FF000080", border="white", add=TRUE)

Tstart.all<-list(D1, D2)

Metrics<-matrix(NA); length(Metrics)<-c(2*1000*4*4); dim(Metrics)<-c(2,1000, 4, 4)
# Metrics: (1) infection pressure, (2) max lice, (3) motile-days, (4) mean lice at 3 sites

for(s in 1:4){ #for each scenario
	
	L<-farmL.all[[s]]
	L<-as.matrix(L)
	dimnames(L)[[2]]<-NULL
	L<-rbind(rep(0, length(T)), L) #add row of zeros for imputing
	
	for(i in 1:1000){
		
		# 1) Total infection pressure
		x1<-seq(-40, 60, dx) #All start at Glendale
		Xind<-match(x1, x)
		
		for(j in 1:2){ #Do two simulations: latter starting 3 weeks earlier
			# T1<-migDate(x1, Xstart=-40, Tstart=Tstart.all[[j]][i], numeric=TRUE)
			# Tind<-findInterval(round(T1, 1)+122, round(T, 1))
			
			# Metrics[j,i,s,1]<-p.lice[1,'phi']*sum(farmL.all[[s]][cbind(Xind, Tind)])*dx	
		
			# # 2) Maximum expected lice
			# Lice_pred<-simulate.lice(p.lice[1,], dist=x1, day=T1+122, log=FALSE)
			
			# Metrics[j,i,s,2]<-max(Lice_pred[[1]][,1]+Lice_pred[[2]][,1]+Lice_pred[[3]][,1])
			
			# # 3) Number of motile-days
			# Metrics[j,i,s,3]<-sum(Lice_pred[[3]][,1])*dt	
			
			# 4) Mean lice at Glacier, Burdwood, Wicklow
			m<-simulate.lice(p.lice[1,], dist=c(32.7, 53, 65), day=round(Tstart.all[[j]][i]+(c(32.7, 53, 65)+40)/p.lice[1,5])+122, log=FALSE)
			Metrics[j,i,s,4]<-mean(m[[1]][,1]+m[[2]][,1]+m[[3]][,1])
			
			# rm(Lice_pred, T1)
			}
				
	}
}

# M1<-list(apply(Metrics[1,,,1], 2, quantile, c(0.5, 0.025, 0.975)), apply(Metrics[2,,,1], 2, quantile, c(0.5, 0.025, 0.975)))
# M2<-list(apply(Metrics[1,,,2], 2, quantile, c(0.5, 0.025, 0.975)), Metrics[2,,,2], 2, quantile, c(0.5, 0.025, 0.975))
# M3<-list(apply(Metrics[1,,,3], 2, quantile, c(0.5, 0.025, 0.975)), Metrics[2,,,3], 2, quantile, c(0.5, 0.025, 0.975))

Metrics.summary<-list(
	InfPressure=list(rbind(apply(Metrics[1,,,1], 2, mean), apply(Metrics[2,,,1], 2, mean)), rbind(apply(Metrics[1,,,1], 2, quantile, c(0.025)), apply(Metrics[2,,,1], 2, quantile, c(0.025))), rbind(apply(Metrics[1,,,1], 2,  quantile, c(0.925)), apply(Metrics[2,,,1], 2,  quantile, c(0.925)))),
	MaxLice=list(rbind(apply(Metrics[1,,,2], 2, mean), apply(Metrics[2,,,2], 2, mean)), rbind(apply(Metrics[1,,,2], 2, quantile, c(0.025)), apply(Metrics[2,,,2], 2, quantile, c(0.025))), rbind(apply(Metrics[1,,,2], 2,  quantile, c(0.925)), apply(Metrics[2,,,2], 2,  quantile, c(0.925)))),
	MotDays=list(rbind(apply(Metrics[1,,,3], 2, mean), apply(Metrics[2,,,3], 2, mean)), rbind(apply(Metrics[1,,,3], 2, quantile, c(0.025)), apply(Metrics[2,,,3], 2, quantile, c(0.025))), rbind(apply(Metrics[1,,,3], 2,  quantile, c(0.925)), apply(Metrics[2,,,3], 2,  quantile, c(0.925)))),
	Monitoring=list(rbind(apply(Metrics[1,,,4], 2, mean), apply(Metrics[2,,,4], 2, mean)), rbind(apply(Metrics[1,,,4], 2, quantile, c(0.025)), apply(Metrics[2,,,4], 2, quantile, c(0.025))), rbind(apply(Metrics[1,,,4], 2,  quantile, c(0.925)), apply(Metrics[2,,,4], 2,  quantile, c(0.925))))
	)

quartz(width=6.3, height=3, pointsize=14)
par(mfrow=c(1,3), mar=c(3,5,2,0), oma=c(1,0,3,1), mgp=c(2.5, 1, 0))

for(i in 1:2){
	bp<-barplot2(Metrics.summary[[i]][[1]], plot.ci=TRUE, ci.l=Metrics.summary[[i]][[2]], ci.u=Metrics.summary[[i]][[3]], las=1, names.arg=LETTERS[1:4], xlab="", ylab=c(expression(integral(L(x,t))), "max [C(x,t)+H(x,t)+M(x,t)]", expression(integral(M(x,t))))[i], beside=TRUE, col=c(4,2))
abline(h=0)
mtext(side=3, line=0.5, c("a) Total infection pressure", "b) Maximum louse load",  "c) Number of motile-days")[i], cex=par('cex'))
if(i==2) legend(-10, max(Metrics.summary[[i]][[3]])*1.4, fill=c(4,2), c("Normal migration timing", "Early migration (e.g., 2010)"), xpd=NA, ncol=2, bty="n")
}

mtext(side=1, outer=TRUE, "Treatment scenario", cex=par('cex'))


##################################################################################################################
# Linking to mortality
##################################################################################################################


c.<-c(0.190, 0.087, 0.299)
Metrics.summary[[5]]<-list(1-exp(-c.[1]*Metrics.summary[[4]][[1]]), 1-exp(-c.[1]*Metrics.summary[[4]][[2]]), 1-exp(-c.[1]*Metrics.summary[[4]][[3]]))

# Predicted mort of fish that migrated out in 2006 (i.e.,returned in 2007):
# 8.3 - 22.3, mean = 15.9 %


quartz(width=6.3, height=2.5, pointsize=14)
par(mfrow=c(2,3), mar=c(3,5,1,0), oma=c(1,0,1,1), mgp=c(2.5, 1, 0), mgp=c(3,1,0))

for(i in c(1:5)){
	bp<-barplot2(Metrics.summary[[i]][[1]], plot.ci=TRUE, ci.l=Metrics.summary[[i]][[2]], ci.u=Metrics.summary[[i]][[3]], las=1, names.arg=LETTERS[1:4], xlab="", ylab=c("Infective larvae", "Number of lice", "Motile-days", "Number of lice", "Percent")[i], beside=TRUE, col=c(grey(0.9), grey(0.3)))
abline(h=0)
mtext(side=3, line=0.5, c("a) Total infection pressure", "b) Maximum louse load",  "c) Number of motile-days", " d) Average louse load", "c) Estimated mortality")[i], cex=par('cex'), adj=0)
	# if(i==5){
		# abline(h=0.159, col=2)
		# abline(h=c(0.083, 0.223), col=2, lty=2)
		# }
	if(i==2)mtext(side=1, "Treatment scenario", cex=par('cex'), line=3)
}
legend(-31.5, 0.18, fill=c(grey(0.9), grey(0.3)), c("Normal", "Early"), bty="n", xpd=NA)
##################################################################################################################
# Lice on fish at Alex's sites
##################################################################################################################
# Use real scenario
L<-farmL.all[[1]]
L<-as.matrix(L)
dimnames(L)[[2]]<-NULL
L<-rbind(rep(0, length(T)), L) #add row of zeros for imputing

loc<-c(32.7, 53, 65)

# Read in Alex's data and bootstrap lice numbers
data<-read.csv("Simulations/Broughton_fishData2006.csv")
data$date<-as.Date(paste(data$year, data$month, data$day, sep="-"))
data$jd<-as.numeric(data$date-as.Date("2006-01-01")+1)
data$chal<-apply(data[,c(12,13)], 1, sum, na.rm=TRUE)
data$mot<-apply(data[,c(14:18)], 1, sum, na.rm=TRUE)
data$cop<-data$Lep_cope; data$cop[is.na(data$cop)]<-0

day1<-rep(unique(data$jd), each=3)
dist1<-rep(c(32.7, 53, 65), length(unique(data$jd)))

# Bootstrap
Lice_obs<-list(); length(Lice_obs)<-9; dim(Lice_obs)<-c(3,3)
for(i in 1:3){for(j in 1:3){Lice_obs[[i,j]]<-matrix(NA, nrow=length(day1), ncol=2)}}

for(i in 1:length(day1)){
	for(j in 1:2){#for species
		dataij<-subset(data, data$location==c("Glacier", "Burdwood", "Wicklow")[as.numeric(as.factor(dist1))[i]]&data$jd==day1[i]&data$species==c("pink", "chum")[j])
		
		for(s in 1:3){#for each stage
			x<-apply(matrix(sample(dataij[,c('cop', 'chal', 'mot')[s]], size=dim(dataij)[1]*1000, replace=TRUE), nrow=1000, ncol=dim(dataij)[1]), 1, mean)
			
			Lice_obs[[s,1]][i,j]<-mean(x, na.rm=TRUE)
			Lice_obs[[s,2]][i,j]<-quantile(x, 0.025, na.rm=TRUE)
			Lice_obs[[s,3]][i,j]<-quantile(x, 0.975, na.rm=TRUE)
			}
		}		
}

Lice_pred<-simulate.lice(p.lice[1,], dist=dist1, day=day1, log=FALSE)

# Plot predcited - observed
par(mfrow=c(3,2), mar=c(3,3,0,0), oma=c(2,2,2,2))
for(s in 1:3){
	for(j in 1:2){
		plotCI(Lice_pred[[1]][,j], Lice_obs[[1,1]][,j], li=Lice_obs[[1,2]][,j], ui=Lice_obs[[1,3]][,j], gap=0, xlab="", ylab="", pch=21, pt.bg="white")
		abline(a=0, b=1, lwd=1.5, lty=2)
		if(s==1) mtext(side=3, c("Pink", "Chum")[j])
	}
	mtext(side=4, c("Cop", "Chal", "Mot")[s], line=0.5)
}
mtext(side=1, outer=TRUE, "Predicted")
mtext(side=2, outer=TRUE, "Observed")

# Doesn't look good. Plot prediction over the whole season at three sites?
day2<-(min(data$jd)-14):(max(data$jd)+14)
Lice_pred2<-list(
		simulate.lice(p.lice[1,], dist=rep(dist1[1], length(day2)), day=day2, log=FALSE),
		simulate.lice(p.lice[1,], dist=rep(dist1[2], length(day2)), day=day2, log=FALSE),
		simulate.lice(p.lice[1,], dist=rep(dist1[3], length(day2)), day=day2, log=FALSE)
	)

par(mfrow=c(3,3), mar=c(2,2,1,1), oma=c(2,2,2,2))
for(i in 1:3){
	for(j in 1:3){
	
		plot(as.Date(paste("2006-", day2, sep=""), format="%Y-%j"), Lice_pred2[[j]][[i]][,1], "l", ylab="", xlab="", ylim=range(c(0, Lice_obs[[i,3]][which(dist1==c(32.7, 53, 65)[j]),1]), na.rm=TRUE))
		lines(as.Date(paste("2006-", day2, sep=""), format="%Y-%j"), Lice_pred2[[j]][[i]][,2], lty=2, lwd=1.5)
		plotCI(unique(data$date), Lice_obs[[i,1]][which(dist1==c(32.7, 53, 65)[j]),1], li=Lice_obs[[i,2]][which(dist1==c(32.7, 53, 65)[j]),1], ui=Lice_obs[[i,3]][which(dist1==c(32.7, 53, 65)[j]),1], gap=0, pch=21, pt.bg="white", add=TRUE)
		plotCI(unique(data$date), Lice_obs[[i,1]][which(dist1==c(32.7, 53, 65)[j]),2], li=Lice_obs[[i,2]][which(dist1==c(32.7, 53, 65)[j]),2], ui=Lice_obs[[i,3]][which(dist1==c(32.7, 53, 65)[j]),2], gap=0, pch=22, pt.bg="#FFFFFF80", add=TRUE)
}}

par(mfrow=c(3,3), mar=c(4,4,2,1))
for(i in 1:3){
	for(j in 1:3){
		ind<-which(dist1==c(32.7, 53, 65)[i])
		ind.data<-which(data$location==c("Glacier", "Burdwood", "Wicklow")[i])
		
		tapply(data[ind.data,c('cop', 'chal', 'mot')[j]], data$date[ind.data], mean)
		
		plot(as.Date(day1[ind], origin="2006-01-01"), Lice_pred[[j]][ind,1], xlab="", ylab="", col=2, las=1)
		points(data$date[ind], data[ind,c('cop', 'chal', 'mot')[j]])
	}
}
