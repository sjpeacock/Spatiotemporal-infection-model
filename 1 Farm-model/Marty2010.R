rm(list=ls())
library(gplots)
setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model")

Z<-read.csv("Farm-model/Marty2010_farm_info.csv")
Z$date<-as.Date(paste("15", Z$month, Z$year, sep="-"), format="%d-%B-%Y")

load("Farm-model/farm_fits_20160209.RData")
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


#############################################################################################
# Number of lice - compare to my data
#############################################################################################

quartz(width=4, height=7)
par(mfrow=c(3,1), mar=c(4,3,1,1), oma=c(2,2,1,1))

for(i in 1:3){
	plot(Z$date[Z$farm_no==c(22,21,11)[i]], Z$motile_Lep[Z$farm_no==c(22,21,11)[i]], xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2, col=2)
	
	if(i<3) plotCI(y.true[[i]][,1], y.true[[i]][,2], li=y.true[[i]][,3], ui=y.true[[i]][,4], add=TRUE, gap=0.3, pch=21, pt.bg="white")
	if(i==3) points(y.true[[i]][,1], y.true[[i]][,2], pch=21, bg="white")
	
	abline(v=treatment.date[i], lty=2)
	abline(v=Z$date[Z$farm_no==c(22,21,11)[i]&is.na(Z$treatment)==FALSE], col=2, lty=2)
	
	}

mtext(side=1, outer=TRUE, "Date (2005/2006)")
mtext(side=2, outer=TRUE, expression(paste("Average motile ",  italic(L.), " ", italic(salmonis),  " per farmed salmon", sep="")))
legend("topright", col=c(1,2), c("My data", "GM data"), pch=1)


#############################################################################################
# Stocking density among farms
#############################################################################################

plot(Z$date[Z$farm_no==22], Z$num_fish[Z$farm_no==22]*10^-3, "l",  xlab="", ylab="", bty="n", ylim=range(Z$num_fish)*10^-3, las=1, cex.axis=1.2)
lines(Z$date[Z$farm_no==21], Z$num_fish[Z$farm_no==21]*10^-3, col=2)
lines(Z$date[Z$farm_no==11], Z$num_fish[Z$farm_no==11]*10^-3, col=3)

# How does this change patterns?
Z$tot_motile_Lep<-Z$motile_Lep*Z$num_fish
quartz(width=4, height=7)
par(mfrow=c(3,1), mar=c(4,4,1,4), oma=c(2,2,1,2))

for(i in 1:3){
	plot(Z$date[Z$farm_no==c(22,21,11)[i]], Z$tot_motile_Lep[Z$farm_no==c(22,21,11)[i]]*10^-3, ylim=range(Z$tot_motile_Lep, na.rm=TRUE)*10^-3, xlab="", ylab="", bty="n", las=1, cex.axis=1.2, col=2)
	
	
	par(new=TRUE)
	if(i<3) plotCI(y.true[[i]][,1], y.true[[i]][,2], li=y.true[[i]][,3], ui=y.true[[i]][,4], gap=0.3, pch=21, pt.bg="white", yaxt="n", xaxt="n", xlim=range(Z$date[Z$farm_no==c(22,21,11)[i]]), ylim=c(0,8), xlab="", ylab="")
	
	if(i==3) plot(y.true[[i]][,1], y.true[[i]][,2], pch=21, bg="white", ylim=c(0,8), yaxt="n", xaxt="n", xlim=range(Z$date[Z$farm_no==c(22,21,11)[i]]), xlab="", ylab="")
	
	axis(side=4, las=1)
	
	abline(v=treatment.date[i], lty=2)
	# abline(v=Z$date[Z$farm_no==c(22,21,11)[i]&is.na(Z$treatment)==FALSE], col=2, lty=2)
	
	}

mtext(side=1, outer=TRUE, "Date (2005/2006)")
mtext(side=2, outer=TRUE, expression(paste("Total motile ",  italic(L.), " ", italic(salmonis),  " on farm", sep="")), col=2)
mtext(side=4, outer=TRUE, expression(paste("Avg. motile ",  italic(L.), " ", italic(salmonis),  " per farmed salmon", sep="")))