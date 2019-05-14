######################################################################################
## Figures for Green's function model paper
######################################################################################
# cols<-c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')
cols<-c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0')
#-------------------------------------------------------------------------------------
## Figure 1: Map
#-------------------------------------------------------------------------------------
setwd("~/Google Drive/Mapping/Broughton Map/")

require(PBSmapping)
load("Data/nepacLLhigh.rda")


data<-read.delim("Data/SeaLice_SummaryDatafile.txt")
names(data)
focal.year<-2006

farms<-read.delim("Data/FarmLocations.txt")
farm.loc<-data.frame(EID=farms$EID,X=farms$lon,Y=farms$lat, name=farms$name, co=farms$company)
farm.loc<-as.EventData(farm.loc, projection="LL")
farm.loc<-farm.loc[c(23,29,30),]
data1<-subset(data, year==focal.year)
X<--data1$lon
Y<-data1$lat
eid<-c(1:length(X))
sampling.loc<-data.frame(EID=eid,X,Y,date=data1$Date, site.name=data1$site.name, temp=data1$temp, sal=data1$sal)
sampling.loc<-sampling.loc[is.na(X)==FALSE,]
sampling.loc<-as.EventData(sampling.loc, projection="LL")

B.box<-data.frame(PID=rep(1,4),POS=c(1:4), X=c(-126.8, -125.5, -125.5, -126.8), Y=c(50.5, 50.5, 51, 51))
B.box<-as.PolySet(B.box, projection="LL")

rivers<-importShapefile("~/Google Drive/Mapping/Broughton Map/Data/bcdrcarto/bcdrcarto_l_1m_v6-0_shp/bcdrcarto_l")
#gpclibPermit()
# Larger Map of BC coast
#quartz(height=4, width=3.5)
pdf("~/Google Drive/Greens/FINALLY/v6/MapBig.pdf", height=4, width=3.5)
plotMap(thinPolys(nepacLLhigh, tol=3), xlim=c(-128.8, -122), ylim=c(48, 53),  col=grey(0.8), bg="white", las=1, xaxt="n", yaxt="n", border=grey(0.4))
#Add box for the Broughton Archipelago
addPolys(B.box, lwd=2, density=0, border=1)

plotMap(nepacLLhigh, xlim=c(-126.8, -125.5), ylim=c(50.5, 51), col=grey(0.8), bg="white", border=grey(0.4), xlab=expression(paste(degree, "Longitude")), ylab=expression(paste(degree, "Latitude")), xaxt="n", yaxt="n")
axis(side=1, at=c(-126.6,-126.4, -126.2, -126.0, -125.8, -125.6), labels=c("-126.6", "-126.4","-126.2", "-126.0", "-125.8","-125.6"), tck=0.03, cex.axis=0.8)
axis(side=1, at=seq(-126.8, -125.5, 0.05), labels=FALSE, tck=0.015)
axis(side=2, at=seq(50.6,50.9,0.1), tck=0.03, las=1, cex.axis=0.8)
axis(side=2, at=seq(50.5, 51, 0.05), labels=FALSE, tck=0.015)

axis(side=3, at=c(-126.6,-126.4, -126.2, -126.0, -125.8, -125.6), labels=FALSE, tck=0.03)
axis(side=3, at=seq(-126.8, -125.5, 0.05), labels=FALSE, tck=0.015)
axis(side=4, at=seq(50.6,50.9,0.1), labels=FALSE, tck=0.03, las=1, cex.axis=0.8)
axis(side=4, at=seq(50.5, 51, 0.05), labels=FALSE, tck=0.015)


#-------------------------------------------------------------------------------------
## Figure 1: Model schematic
#-------------------------------------------------------------------------------------
# In LaTeX via tikz

#-------------------------------------------------------------------------------------
## Figure 3: The average number of motile L. salmonis per farmed salmon 
## on three salmon farms under four different treatment scenarios:
#-------------------------------------------------------------------------------------
load("Workspace/Scenarios_20181130.RData")
pdf(file="Figures/Fig3.pdf", width=5.3, height=4, pointsize=10)

par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(0,1,0,0))
topText<-c("a) Scenario A", "b) Scenario B", "c) Scenario C", "d) Scenario D")
for(i in 1:4){
	plot(as.Date(times, origin="2005-09-01"), predicted.f[[i]][,1], "l",col=1, ylim=c(0, 8), bty="l", las=1, xlab="", ylab="")#max(predicted.f[[i]])
	polygon(x=c(as.Date("2006-04-01"), as.Date("2006-04-01"), as.Date("2006-06-30"), as.Date("2006-06-30")), y=c(-1, 10, 10, -1), border=NA, col="#00000020")
	lines(as.Date(times, origin="2005-09-01"), predicted.f[[i]][,2], lty=2)
	lines(as.Date(times, origin="2005-09-01"), predicted.f[[i]][,3], lty=3)
	#abline(v=as.Date(tT[[i]], origin="2005-09-01"), lty=c(1:3))
	abline(h=3, col=cols[1])#, col=grey(0.8), lwd=2)
	mtext(side=3, adj=0, line=0.5, topText[i])
	if(i==2) legend("topleft", lty=c(1,2,3), c("Farm 1", "Farm 2", "Farm 3"), bg="white", bty="n")#legend(as.Date("2005-09-01"), 13, lwd=c(1,1,1,1,15), col=c(rep(1,3), 2,grey(0.8)), lty=c(1,2,3,1,1), c("Sargeaunts", "Humphrey", "Burdwood", "Treatment threshold", "Juvenile salmon migration"), bg="white", ncol=3, xpd=NA, bty="n")
}
mtext(side=2, outer=TRUE, expression(paste("Average number of motile ", italic(L), ". ",  italic(salmonis), " per farmed salmon")), line=-0.5)

dev.off()

#-------------------------------------------------------------------------------------
## Figure 4: Farm predictions
#-------------------------------------------------------------------------------------
require(gplots)
load("Workspaces/farmModels_13June2012.RData")
load("Workspaces/farm_fits_20160209.RData")

pdf(file="Figures/Fig4.pdf", width=4, height=7)
par(mfrow=c(3,1), mar=c(4,3,1,1), oma=c(2,2,1,1))

plot(as.Date(T.all, origin="1970-01-01"), fset(p.fitted[[1]][,1], T.all, T.slice.all[1]), "l", xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2)
polygon(x=c(T.all, rev(T.all)), y=c(CI.pred[[1]][,1], rev(CI.pred[[1]][,2])), border=NA, col="#00000030")
plotCI(as.Date(unique(Z$Date[Z$Farm=="SP"]), origin="1970-01-01"), Lice.summary[[1]][,1], ui=Lice.summary[[1]][,3], li=Lice.summary[[1]][,2], add=TRUE)
abline(v=T.slice.all[1], lty=2)
mtext(side=3, adj=0, "a) Farm 1", line=0.5)

plot(as.Date(T.all, origin="1970-01-01"), fset(p.fitted[[2]][,1], T.all, T.slice.all[2]), "l", xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2)
polygon(x=c(T.all, rev(T.all)), y=c(CI.pred[[2]][,1], rev(CI.pred[[2]][,2])), border=NA, col="#00000030")
plotCI(as.Date(unique(Z$Date[Z$Farm=="HR"]), origin="1970-01-01"), Lice.summary[[2]][,1], ui=Lice.summary[[2]][,3], li=Lice.summary[[2]][,2], add=TRUE)
abline(v=T.slice.all[2], lty=2)
mtext(side=3, adj=0, "b) Farm 2", line=0.5)

plot(as.Date(T.all, origin="1970-01-01"), fset(p.fitted[[3]][,1], T.all, T.slice.all[3]), "l", xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2)
polygon(x=c(T.all, rev(T.all)), y=c(CI.pred[[3]][,1], rev(CI.pred[[3]][,2])), border=NA, col="#00000030")
points(as.Date(Z$Date[Z$Farm=="BG"], origin="1970-01-01"), Z$Lice[Z$Farm=="BG"]/20)
abline(v=T.slice.all[3], lty=2)
mtext(side=3, adj=0, "c) Farm 3", line=0.5)

mtext(side=2, outer=TRUE, "Average lice per farm salmon")
mtext(side=1, outer=TRUE, "Date (2005/2006)")
dev.off()


#-------------------------------------------------------------------------------------
## Figure 5: The simulated densities of infectious copepodites
#-------------------------------------------------------------------------------------
load("Workspaces/Scenarios_20181130.RData")

source("filledContour.R")

pal<-colorRampPalette(c("white", 1))
source("4 Simulations/sim-grid.R")

pdf(file="Figures/SimDens.pdf", width=6.3, height=6, pointsize=10)

par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(2,2,1,0))
for(i in 1:4){
	filledContour(farmL.all[[i]][seq(1, 280, 2), seq(1, 606, 4)], x=x[seq(1, 280, 2)], y=as.Date(T, origin="2005-09-01")[seq(1, 606, 4)], zlim=c(0,0.05), color.palette=pal, xlab="", ylab="")
	contour(farmL.all[[i]], x=x, y=as.Date(T,origin="2005-09-01"), add=TRUE)
	mtext(side=3, line=0.05, paste("Scenario ", LETTERS[i]))
	arrows(x0=-40, y0=as.Date(paste("2006-", e1, sep=""), format="%Y-%j"), x1=60, y1=migDate(60, Xstart=-40, Tstart=e1), col=cols[5], length=0.08, xpd=NA, lwd=2)
	points(-40, as.Date(paste("2006-", e1, sep=""), format="%Y-%j"), pch=19, col=cols[5])
	arrows(x0=-40, y0=as.Date(paste("2006-", e2, sep=""), format="%Y-%j"), x1=60, y1=migDate(60, Xstart=-40, Tstart=e2), col=cols[1], length=0.08, xpd=NA)
	points(-40, as.Date(paste("2006-", e2, sep=""), format="%Y-%j"), pch=21, col=cols[1], bg="white")
	
	
	if(i>=3) abline(h=as.Date(tT[[i]][1], origin="2005-09-01"), lty=2)
	abline(v=c(-3.7, 4.0, 53), lty=c(1:3), lwd=1.2)	
	
}
mtext(side=1, outer=TRUE, "Distance along migration (km)")
mtext(side=2, outer=TRUE, "Date", las=0)

dev.off()
#-------------------------------------------------------------------------------------
## Figure 6: Simulation metrics
#-------------------------------------------------------------------------------------
load("Workspaces/Scenarios_20181130.RData")

pdf(file="Figures/Fig6.pdf", width=6.3, height=2.5)
par(mfrow=c(1,3), mar=c(3,5,1,0), oma=c(1,0,1,1), mgp=c(2.5, 1, 0), mgp=c(3,1,0))

for(i in c(c(1,2,5))){
	bp<-barplot2(Metrics.summary[[i]][[1]], plot.ci=TRUE, ci.l=Metrics.summary[[i]][[2]], ci.u=Metrics.summary[[i]][[3]], las=1, names.arg=LETTERS[1:4], xlab="", ylab=c("Total infection pressure\nalong migration route", "Max. sea lice per fish", "Motile-days", "Number of lice", "Estimated mortality of wild salmon")[i], beside=TRUE, col=cols[c(4,1)], ci.width=0.4)
	abline(h=0)
	mtext(side=3, line=0.5, c("a)", "b)",  "c) Number of motile-days", " d) Average louse load", "c) ")[i], cex=par('cex'), adj=0)
	# if(i==5){
	# abline(h=0.159, col=2)
	# abline(h=c(0.083, 0.223), col=2, lty=2)
	# }
	if(i==2)mtext(side=1, "Treatment scenario", cex=par('cex'), line=3)
	if(i==1)legend(5.5, 4.1, fill=cols[c(4,1)], c("Normal", "Early"), bty="n", xpd=NA, title="Migration timing")
}

dev.off()

######################################################################################
# SUPPORTING INFROMATION #
######################################################################################

#-------------------------------------------------------------------------------------
## SI Figure 1: Farm predictions
#-------------------------------------------------------------------------------------
load("Code/Cope distribution/CopedistWorkspace_June14.RData")
p<-read.csv("Code/JAGS/2006/neg binom/dcloneFittedParameters.csv")

t.sim0<-seq(50,150,2)
t.ind<-match(t.sim0, T);t.ind[length(t.ind)]<-2800

ani.options(outdir = "~/Documents/Greens/FINALLY/Figures/FarmVambient/3SeparateFarms")
saveLatex({

## put any code here to produce several plots
	for(i in 1:length(t.sim0)){
		plot(x, p$k[1]*p$phi[1]*copedist[,t.ind[i]], "l", bty="l", las=1, ylab="Planktonic copepidites", xlab="Distance", ylim=c(0, max(p$k[1]*p$phi[1]*copedist)), lwd=2, col=grey(0.8))
		lines(x, p$k[1]*p$phi[1]*SP[,t.ind[i]]/scale.par)
		lines(x, p$k[1]*p$phi[1]*HR[,t.ind[i]]/scale.par)
		lines(x, p$k[1]*p$phi[1]*BW[,t.ind[i]]/scale.par)
		lines(x, rep(p$k[1], length(x)), lty=2)
		abline(v=c(-3.71, 4, 53), lty=3)
		mtext(side=3, as.Date("2006-01-01")+t.sim0[i])
		#if(t.sim0[i]>=67) lines(x, p$k[1]*p$phi[1]*copedist[,which(T==67)], col="#00000040", lwd=1.2)
		}
	
	},
interval = 0.4, ani.dev = 'pdf', ani.type = 'pdf', ani.height = 4, ani.width = 5, ani.opts='controls,width=4in', pdflatex = '/usr/texbin/pdflatex', overwrite=TRUE, documentclass = paste("\\documentclass{article}", "\\usepackage[margin=0.3in]{geometry}", sep="\n"))

#-------------------------------------------------------------------------------------
## Observed and predicted surfaces
#-------------------------------------------------------------------------------------
library(here)
library(dclone)
library(parallel)
library(boot) # for inv.logit function

# Load in results from model fitting by "fit-model-JAGS.R"
load('Workspaces/LiceFits2Data_20160308.RData')

# Summaryize MCMC output
S<-summary(mod)

source('3 Fitting/liceBoot.R')
source("3 Fitting/sim-model.R")
d<-4 # Scale of grid for plotting
n.dist<-length(seq(-44,68,d))
n.day<-length(seq(100,144,d))
dist2<-rep(seq(-44,68,d), n.day)
day2<-rep(seq(100,144,d), each=n.dist)
L<-farmL

sim.out<-simulate.lice(p=S[[1]][,1], dist=dist2, day=day2)

j<-1 #species = pink

# colour scheme
colPred <- c(surf = "#FFFFFF80", grid = "#d7191c80")
# colDat <- c(above = "#2c7bb6", below = "#abd9e9")
colDat <- c(above = 1, below = "#00000030")
# pdf(file = "Figures/Supplement-PinkFits2Data.pdf", width=4, height=9)


colPred <- c(surf = "#FFFFFF80", grid = grey(0.6))
colDat <- c(above = "#2c7bb6", below = "#d7191c60")

# quartz(width = 5, height = 5, pointsize=10)
# par(mfrow=c(1,1), mar=c(1,0,1,0), oma=c(0,0,1,0))

quartz(width = 4, height = 9, pointsize=10)
par(mfrow=c(3,1), mar=c(1,0,1,0), oma=c(0,0,1,0))

for(i in 1:3){
	Z<-matrix(sim.out[[i]][,j], nrow=n.dist, ncol=n.day, byrow=FALSE)
	pmat1<-persp(z=Z, x=seq(-44,68,d), y=seq(100,144,d), xlab="Distance", ylab="Time", zlab=c("C(x,t)", "H(x,t)", "M(x,t)")[i], theta=-140, phi=15, col=colPred['surf'], zlim=c(0, max(Lice.boot[i,,j,1])), border = colPred['grid'], ticktype = "detailed")
	
	Lmean<-trans3d(x=dist, y=day, z=Lice.boot[i,,j,1], pmat=pmat1)
	Lmin<-trans3d(x=dist, y=day,  z=Lice.boot[i,,j,2], pmat=pmat1)
	Lmax<-trans3d(x=dist, y=day,  z=Lice.boot[i,,j,3], pmat=pmat1)
	
	for(k in 1:length(day)){
		ind<-c(findInterval(dist[k],seq(-44,68,d)), findInterval(day[k],seq(100,144,d)))
		
		# text(Lmean$x[k], Lmean$y[k], k, xpd=NA)
		
		if(Lice.boot[i,k,j,1]>=Z[ind[1], ind[2]]){ #if obs>pred
			points(Lmean$x[k], Lmean$y[k], pch=19, xpd=NA, col=colDat['above'])
			segments(Lmean$x[k], Lmean$y[k], Lmax$x[k], Lmax$y[k], col=colDat['above'])
			
			if(Lice.boot[i,k,j,2]>Z[ind[1], ind[2]]){
				segments(Lmean$x[k], Lmean$y[k], Lmin$x[k], Lmin$y[k], col=colDat['above'])
			}else{
				modMin<-trans3d(x=dist[k], y=day[k], z=Z[ind[1], ind[2]], pmat=pmat1)
				segments(Lmean$x[k], Lmean$y[k], modMin$x, modMin$y, col=colDat['above'])
				segments(modMin$x, modMin$y, Lmin$x[k], Lmin$y[k], col= colDat['below'])
			}
		} # end if obs > pred
		
		if(Lice.boot[i,k,j,1]<Z[ind[1], ind[2]]){ #if obs < pred	
			points(Lmean$x[k], Lmean$y[k], pch=19, col=colDat['below'])
			segments(Lmean$x[k], Lmean$y[k], Lmin$x[k], Lmin$y[k], col=colDat['below'])
			
			if(Lice.boot[i,k,j,3]<Z[ind[1], ind[2]]){
				segments(Lmean$x[k], Lmean$y[k], Lmax$x[k], Lmax$y[k], col=colDat['below'])
			}else{
				modMin<-trans3d(x=dist[k], y=day[k], z=Z[ind[1], ind[2]], pmat=pmat1)
				segments(Lmean$x[k], Lmean$y[k], modMin$x, modMin$y, col=colDat['below'])
				segments(modMin$x, modMin$y, Lmax$x[k], Lmax$y[k], col=colDat['above'])
			}
		} #end if obs < pred
		
	} # end k
	
	# par(new=TRUE)
	# persp(z=Z, x=seq(-44,68,d), y=seq(100,144,d), xlab="", ylab="", zlab="", theta=-140, phi=15, col=colPred['surf'], zlim=c(0, max(Lice.boot[i,,j,1])), border = colPred['grid'])
	
	mtext(side=3, adj=0, c("a) Copepodid", "b) Chalimus", "c) Motile")[i])	
	
} #end stage i

#-------------------------------------------------------------------------------------
## Observed versus predicted over time
#-------------------------------------------------------------------------------------
colR <- colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))

sim.out2<-simulate.lice(p=S[[1]][,1], dist=dist, day=day)

for(i in 1:3){
	plotCI(sim.out2[[i]][,j], Lice.boot[i,,j,1], li = Lice.boot[i,,j,2], ui = Lice.boot[i,,j,3], gap=0, pch=21, pt.bg = colR(n=length(unique(day)))[match(1:length(unique(day)), order(day))], xlim=c(0.12, 0.15), ylim=c(0, 0.3))
	abline(a =0, b=1, lty=2)
#-------------------------------------------------------------------------------------
## Temperature and salinity over course of sampling
#-------------------------------------------------------------------------------------

siteDat <- read.csv("Figures/site_data2006.csv")

index <- read.csv("3 Fitting/index.csv")

head(siteDat)

range(siteDat$temp, na.rm=TRUE)
range(siteDat$sal, na.rm=TRUE)

mean(siteDat$temp, na.rm=TRUE)

# How many sites sampled in a day?
range(tapply(siteDat$temp, siteDat$date, length))

siteDat$date <- as.Date(siteDat$FullDate, format = "%Y-%m-%d")

plot(siteDat$date, siteDat$temp)

meanTemp <- tapply(siteDat$temp, siteDat$date, mean)
sdTemp <- tapply(siteDat$temp, siteDat$date, sd)
meanSal <- tapply(siteDat$sal, siteDat$date, mean)
sdSal <- tapply(siteDat$sal, siteDat$date, sd)

tempRange <- seq(7, 15, 0.1)
devTimeC <- tau(tempRange, beta[1,1], beta[1,2])
tPoints <- rev(tempRange)[findInterval(c(2:5), rev(devTimeC))]
plot(tempRange, devTimeC, "l")

quartz(width = 4.5, height = 6, pointsize =10)

par(mfrow=c(2, 1), mar=c(3,4,2,4.5), oma=c(2,0,0,0))
plotCI(as.Date(names(meanTemp), format = "%Y-%m-%d"), meanTemp, type = "n", xlab="", ylab=expression(paste("Surface temperature (", degree, "C)")), las=1, liw=sdTemp, uiw=sdTemp, lwd = NA, bty="u")
axis(side = 4, at = tPoints, labels = c(2:5), las=1)
mtext(side = 4, "Development time\nof copepodites (days)", line=3)
abline(h = 10, lty=2, col = "#d7191c")
abline(h = mean(siteDat$temp, na.rm=TRUE), lwd = 2, col = "#2c7bb6")
plotCI(as.Date(names(meanTemp), format = "%Y-%m-%d"), meanTemp, liw=sdTemp, uiw=sdTemp, gap=0, pch=21, pt.bg="white", sfrac = 0.008, add=TRUE)
mtext(side = 3, adj=0, line = 0.5, "a)")
legend("topleft", lwd=c(2, 1), lty=c(1,2), col=c("#2c7bb6","#d7191c"), c("mean", "assumed value"), bty="n")

plotCI(as.Date(names(meanTemp), format = "%Y-%m-%d"), meanSal, type ="n", xlab="", ylab="Surface salinity (ppm)", las=1, liw=sdSal, uiw=sdSal, lwd = NA, bty="u")
# axis(side = 4, at = sPoints, labels = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8), las=1, xpd=NA)
# mtext(side = 4, "Mortality rate of copepodites", line=3)
abline(h = 30, lty=2, col = "#d7191c")
abline(h = mean(siteDat$sal, na.rm=TRUE), lwd = 2, col = "#2c7bb6")
plotCI(as.Date(names(meanTemp), format = "%Y-%m-%d"), meanSal, liw=sdSal, uiw=sdSal, gap=0, pch=21, pt.bg="white", sfrac = 0.008, add=TRUE)
mtext(side = 3, adj=0, line = 0.5, "b)")

# What impact on parameters
#Fixed parameters from Stien et al. (2005)
tauC <- 0.525^-2 # Development time of copepodids
tauH <- 0.250^-2 # Development time of chalimus
tauM <- 0.187^-2 # Development time of motiles

beta <- matrix(c(24.79, 74.70, 67.47, 0.525, 0.250, 0.187), nrow = 3, ncol = 2, dimnames = list(c("C", "H", "M")))

tau <- function(temp, beta1, beta2){
	(beta1/(temp - 10 + beta1*beta2))^2
}

tau(c(7.7, 12.5), beta[1,1], beta[1,2])
tau(c(7.7, 12.5), beta[2,1], beta[2,2])
tau(c(7.7, 12.5), beta[3,1], beta[3,2])


tau(10, beta[1,1]+c(-1.96, 1.96)*1.43, beta[1,2]+c(-1.96, 1.96)*0.017)

# Salinity relationship from Groner et al. (2016)
# time (hrs) = exp(beta0 + beta1*salinity)
M <- function(sal, denom, exponent){
	1/(1+(sal/denom)^exponent)
}


M(mean(siteDat$sal, na.rm=TRUE), 19.09, 7.11)

salRange <- seq(10, 40, 0.1)
Ms <- M(salRange, 19.09, 7.11)
sPoints <- rev(salRange)[findInterval(c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8), rev(Ms))]

###################

quartz(width = 5, height = 6, pointsize = 12)
layout(matrix(c(1,1,2,3,3,4), nrow=2, byrow=TRUE))
par(oma=c(2,0,0,0))

par(mar=c(4,4,2,0))
plotCI(as.Date(names(meanTemp), format = "%Y-%m-%d"), meanTemp, type = "n", xlab="", ylab=expression(paste("Surface temperature (", degree, "C)")), las=1, liw=sdTemp, uiw=sdTemp, lwd = NA)
abline(h = 10, lty=2, col = "#d7191c")
abline(h = mean(siteDat$temp, na.rm=TRUE), lwd = 2, col = "#2c7bb6")
abline(h = range(meanTemp, na.rm=TRUE), lty = 2, col = "#2c7bb6")
plotCI(as.Date(names(meanTemp), format = "%Y-%m-%d"), meanTemp, liw=sdTemp, uiw=sdTemp, gap=0, pch=21, pt.bg="white", sfrac = 0.008, add=TRUE)
mtext(side = 3, adj=0, line = 0.5, "a)")
legend("topleft", lwd=c(2, 1,1), lty=c(1,2,2), col=c("#2c7bb6","#2c7bb6","#d7191c"), c("overall mean", "range in daily avg.", "assumed constant value"), bty="n")

par(mar=c(4,0,2,4))
u <- par('usr')
plot(tau(seq(u[3], u[4], 0.01), beta[1,1], beta[1,2]), seq(u[3], u[4], 0.01), "l", lwd=2, yaxs="i", yaxt="n", xlab="Development time\n of copepodids (days)")
u1 <- par('usr')
segments(x0 = u1[1], x1=tauC, y0 = 10, y1 = 10, lty=2, col = "#d7191c")
segments(x0 = tauC, x1=tauC, y0 = 10, y1 = u1[3], lty=2, col = "#d7191c")

segments(x0 = rep(u1[1],2), x1=tau(range(meanTemp, na.rm=TRUE), beta[1,1], beta[1,2]), y0 = range(meanTemp, na.rm=TRUE), y1 = range(meanTemp, na.rm=TRUE), lty=2, col = "#2c7bb6")
segments(x0 = tau(range(meanTemp, na.rm=TRUE), beta[1,1], beta[1,2]), x1=tau(range(meanTemp, na.rm=TRUE), beta[1,1], beta[1,2]), y0 = range(meanTemp, na.rm=TRUE), y1 = u1[3], lty=2, col = "#2c7bb6")
text(4.4, 15, "Stien et al. (2005)")
axis(side=4, las=1)

par(mar=c(4,4,2,0))
plotCI(as.Date(names(meanSal), format = "%Y-%m-%d"), meanSal, type ="n", xlab="", ylab="Surface salinity (ppm)", las=1, liw=sdSal, uiw=sdSal, lwd = NA)
# axis(side = 4, at = sPoints, labels = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8), las=1, xpd=NA)
# mtext(side = 4, "Mortality rate of copepodites", line=3)
abline(h = 30, lty=2, col = "#d7191c")
abline(h = mean(siteDat$sal, na.rm=TRUE), lwd = 2, col = "#2c7bb6")
abline(h = range(meanSal[meanSal>15], na.rm=TRUE), lty = 2, col = "#2c7bb6")
plotCI(as.Date(names(meanTemp), format = "%Y-%m-%d"), meanSal, liw=sdSal, uiw=sdSal, gap=0, pch=21, pt.bg="white", sfrac = 0.008, add=TRUE)
mtext(side = 3, adj=0, line = 0.5, "b)")
points(as.Date(names(meanSal[meanSal < 15]), format = "%Y-%m-%d"), meanSal[meanSal < 15], pch=19)

par(mar=c(4,0,2,4))
u <- par('usr')
S <- seq(u[3], u[4], 0.1)
plot(1 - 1/(1+(S/21.22)^5.82), S, "l", lwd=2, yaxs="i", yaxt="n", xlab="Proportion of\ncopepodids attaching")
u1 <- par('usr')

Srange <- range(meanSal[meanSal>15], na.rm=TRUE)
segments(x0 = rep(u1[1],2), x1=1 - 1/(1+(Srange/21.22)^5.82), y0 = Srange, y1 = Srange, lty=2, col = "#2c7bb6")
segments(x0 = 1 - 1/(1+(Srange/21.22)^5.82), x1=1 - 1/(1+(Srange/21.22)^5.82), y0 = Srange, y1 = u1[3], lty=2, col = "#2c7bb6")

text(0.42, 33.5, "Groner et al. (2016)")
