rm(list=ls())

library(dclone); library(gplots)

setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model")
load('0 Original model/OriginalModel_DClone_20161109.RData')

S<-summary(mod)

#######################################################
# table of estimates for Sweave
#######################################################
p.lice.log<-matrix(c(S[[1]][,1], S[[1]][,1]-1.96*S[[1]][,3], S[[1]][,1]+1.96*S[[1]][,3]), nrow=3, ncol=8, byrow=TRUE)


p.lice<-cbind(plogis(p.lice.log[,1:4]), exp(p.lice.log[,5:8]))
rownames(p.lice)<-c("mean", "2.5%", "97.5%")
colnames(p.lice)<-c("Csc", "Csh", "Psc", "Psh", "V",   "k" ,  "phi", "r")


# #######################################################
# # Compare to priors
# #######################################################
# prior<-data.frame(par.names, mu, sig)

# par(mfrow=c(2,4), mar=c(4,4,2,1))
# for(i in 1:8){
	# xR<-range(c(mod[[1]][,i], mod[[2]][,i], mod[[3]][,i], init.p[[1]][which(names(init.p[[1]])==colnames(mod[[1]])[i])], init.p[[2]][which(names(init.p[[2]])==colnames(mod[[1]])[i])], init.p[[3]][which(names(init.p[[3]])==colnames(mod[[1]])[i])]))
	# hist(c(mod[[1]][,i], mod[[2]][,i], mod[[3]][,i]), border=NA, xlab="parameter value", ylab="density", freq=FALSE, main=colnames(mod[[1]])[i], xlim=xR)
	# for(j in 1:3) lines(density(mod[[j]][,i]), col=j, lwd=2)
	# xP<-seq(xR[1], xR[2], length.out=100)
	# lines(xP, dnorm(xP, mean=prior$mu[which(prior$par.names==colnames(mod[[1]])[i])], sd=prior$sig[which(prior$par.names==colnames(mod[[1]])[i])]), col=4)
	# for(j in 1:3) abline(v=init.p[[j]][which(names(init.p[[j]])==colnames(mod[[1]])[i])], col=j, lty=2)
# }

# initial values don't seem to have much to do with it. Maybe just run waaaay longer?

#######################################################
# Data cloning results
#######################################################
DC<-dctable(mod)
k<-c(1:3,5:8)


quartz(width=6.3, height=3, pointsize=11)
par(mfrow=c(2,4), mar=c(1,1,0,0), oma=c(4,4,1,1))
for(i in 1:15){
	plot(DC[[i]]$n.clones[k], DC[[i]]$sd[k]^2/(DC[[i]]$sd[1]^2), "l", bty="l", las=1, ylab="", xlab="", xaxt="n", yaxt="n", ylim=c(0,1))
	lines(1:10, 1/(1:10), lty=2)
	points(DC[[i]]$n.clones[k], DC[[i]]$sd[k]^2/(DC[[i]]$sd[1]^2), pch=21, bg="white")
	mtext(side=3, colnames(mod[[1]])[i], line=-2)
	if(i>4) axis(side=1)else axis(side=1, labels=FALSE)
	if(i==1|i==5) axis(side=2, las=1) else axis(side=2, labels=FALSE)
}
mtext(side=1, outer=TRUE, "Number of clones", line=2)
mtext(side=2, outer=TRUE, "Scaled variance", line=2)


par(mfrow=c(2,3), mar=c(4,4,2,1))
for(i in 6:7){
	for(j in 8:10){
		plot(as.numeric(c(mod[[1]][,i], mod[[2]][,i], mod[[3]][,i])), as.numeric(c(mod[[1]][,j], mod[[2]][,j], mod[[3]][,j])), xlab=colnames(mod[[1]])[i], ylab=colnames(mod[[1]])[j], bty="l", col="#00000030")}}


######################################################################################
#
# For presentation: lambda_c (i=2,3), lh_R (4), Lm_R (5), sc and sh (12-15)
######################################################################################
######################################################################################

par.names<-c(expression(italic(D)), expression(paste(lambda[c], "-pink")), expression(paste(lambda[c], "-chum")), expression(lambda[h]/lambda[c]), expression(lambda[m]/lambda[c]), expression(paste(italic(k), "-pink")), expression(paste(italic(k), "-chum")), expression(phi[1]), expression(phi[2]), expression(phi[3]), expression(italic(r)), expression(paste(italic(s[c]), "-pink")), expression(paste(italic(s[c]), "-chum")), expression(paste(italic(s[h]), "-pink")), expression(paste(italic(s[h]), "-chum")))


quartz(width=5, height=4)

par(mfrow=c(3,2), mar=c(2,2,2,1), oma=c(2,2,0,0))
for(i in c(4:5,12:15)){
	y<-(DC[[i]]$sd)^2/(DC[[i]]$sd[1])^2
	plot(DC[[i]]$n.clones, y, "o", pch=21, bg="white", ylim=c(0, 1), las=1, bty="n", xlim=c(0,11), xlab="", yaxt="n", ylab="")
	axis(side=2, at=c(0, 0.5, 1), las=1)
	lines(DC[[i]]$n.clones, 1/DC[[i]]$n.clones, lty=3)
	points(DC[[i]]$n.clones[which(DC[[i]]$r.hat>1.1)], y[which(DC[[i]]$r.hat>1.1)], pch=19, col=2, cex=0.75)
	mtext(side=3, par.names[i])
	}
mtext(side=1, outer=TRUE, "Number of clones", cex=par('cex'), line=1)
mtext(side=2, outer=TRUE, "Scaled variance", cex=par('cex'), line=1)


quartz(width=3, height=3, pointsize=10)
plot(as.numeric(c(mod[[1]][,5], mod[[2]][,5], mod[[3]][,5])), as.numeric(c(mod[[1]][,14], mod[[2]][,14], mod[[3]][,14])), bty="n", xlab=expression(paste("logit ", lambda[m]/lambda[c])), ylab=expression(paste("log ", s[h], "-pink")), col="#00000030", las=1)


quartz(width=3, height=4, pointsize=10)
par(mfrow=c(2,1), mar=c(4,2,2,1), oma=c(0,0,0,0))
hist(inv.logit(c(mod[[1]][,12], mod[[2]][,12], mod[[3]][,12])), xlab=expression(paste(s[c], "-pink")), freq=FALSE, yaxt="n", col=2, border="white", main="")
arrows(par('usr')[1], par('usr')[3], par('usr')[1], par('usr')[4], length=0.08, xpd=NA)
mtext(side=2, "Density", line=1)
hist(inv.logit(c(mod[[1]][,14], mod[[2]][,14], mod[[3]][,14])), xlab=expression(paste(s[h], "-pink")), freq=FALSE, , yaxt="n", col=2, border="white", main="")
arrows(par('usr')[1], par('usr')[3], par('usr')[1], par('usr')[4], length=0.08, xpd=NA)
mtext(side=2, "Density", line=1)

y<-exp(c(mod[[1]][,5], mod[[2]][,5], mod[[3]][,5]))
hist(y, xlab=par.names[5], freq=FALSE, , yaxt="n", col=2, border="white", main="", ylab="")
lines(dlnorm(seq(min(y), max(y), length.out=100), meanlog=mu[8], sig[8]), lty=3)
arrows(par('usr')[1], par('usr')[3], par('usr')[1], par('usr')[4], length=0.08, xpd=NA)
mtext(side=2, "Density", line=1)

######################################################################################
######################################################################################

# Survival correlations
par(mfrow=c(15, 15), mar=c(0,0,0,0), oma=c(0,2,2,0))
for(i in 1:15){
	for(j in 1:15){
		plot(as.numeric(c(mod[[1]][seq(1,2000,10),i], mod[[2]][seq(1,2000,10),i], mod[[3]][seq(1,2000,10),i])), as.numeric(c(mod[[1]][seq(1,2000,10),j], mod[[2]][seq(1,2000,10),j], mod[[3]][seq(1,2000,10),j])), xlab="", ylab="", bty="o", col="#00000030", xaxt="n", yaxt="n")
		if(i==1) mtext(side=3, colnames(mod[[1]])[j], cex=0.6)
		if(j==1) mtext(side=2, colnames(mod[[1]])[i], cex=0.6)
		
		}}

# Correlation betwen sh and Lm_R
		
		
		
######################################################################################
# Parameters of interest
######################################################################################

est<-exp(rbind(HPDinterval(as.mcmc(c(mod[[1]][,6], mod[[2]][,6], mod[[2]][,6]))), HPDinterval(as.mcmc(c(mod[[1]][,7], mod[[2]][,7], mod[[2]][,7]))), 
HPDinterval(as.mcmc(c(mod[[1]][,8], mod[[2]][,8], mod[[2]][,8]))),
HPDinterval(as.mcmc(c(mod[[1]][,9], mod[[2]][,9], mod[[2]][,9]))),
HPDinterval(as.mcmc(c(mod[[1]][,10], mod[[2]][,10], mod[[2]][,10])))))
rownames(est)<-c("kp", "kc", "phi1", "phi2", "phi3")



####################################################################################################
## Look at MCMC output
####################################################################################################
thin<-1
layout(matrix(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10,11,11,12,13,13,14), nrow=7, ncol=3, byrow=TRUE))
par(mar=c(2,1,0,1), oma=c(2,4,3,0), mgp=c(3.5, 1, 0))
for(i in 1:7){
	plot(mod[[1]][seq(1, dim(mod[[1]])[1], thin),i], type="l", lty=2, col=1, ylab=par.names[i], xlab="", las=1, xpd=NA, bty="l", xaxt="n")
	if(i==7){axis(side=1); mtext(side=1, "Iteration", line=2.5, cex=0.8)} else axis(side=1, labels=FALSE)
	lines(mod[[2]][seq(1, dim(mod[[1]])[1], thin),i], lty=2, col=2)
	lines(mod[[3]][seq(1, dim(mod[[1]])[1], thin),i], lty=2, col=3)
	if(i==1) mtext(side=3, line=1, "MCMC trace", cex=0.8)
	
	plot(density(mod[[1]][,i]), bty="n", main="", ylab="", las=1, yaxt="n")
	arrows(par()$usr[1],par()$usr[3], par()$usr[1], par()$usr[4], length=0.08, xpd=NA) 
	lines(density(mod[[2]][,i]), col=2); lines(density(mod[[3]][,i]), col=3, ylab="", xlab=colnames(mod[[1]])[i])
	if(i==7) mtext(side=1, "Parameter value", line=2.5, cex=0.8)
	if(i==1) mtext(side=3, line=1, "Density", cex=0.8)
	}
}

#######################################################
# Plot model results
#######################################################
source('Fitting/liceBoot.R')
source("Fitting/sim-model.R")
d<-4 # Scale of grid for plotting
n.dist<-length(seq(-44,68,d))
n.day<-length(seq(100,144,d))
dist2<-rep(seq(-44,68,d), n.day)
day2<-rep(seq(100,144,d), each=n.dist)
L<-farmL

sim.out<-simulate.lice(p=S[[1]][,1], dist=dist2, day=day2)

j<-1 #species
i<-1 #stage 
quartz(width=4, height=9)
par(mfrow=c(3,1), mar=c(1,0,1,0), oma=c(0,0,1,0))
for(i in 1:3){
	Z<-matrix(sim.out[[i]][,j], nrow=n.dist, ncol=n.day, byrow=FALSE)
	pmat1<-persp(z=Z, x=seq(-44,68,d), y=seq(100,144,d), xlab="Distance", ylab="Time", zlab=c("C(x,t)", "H(x,t)", "M(x,t)")[i], theta=-150, phi=25, col="#FFFFFF30", zlim=c(0, max(Lice.boot[i,,j,1])))#, main=c("Pink", "Chum")[j]
	
	Lmean<-trans3d(x=dist, y=day, z=Lice.boot[i,,j,1], pmat=pmat1)
	Lmin<-trans3d(x=dist, y=day,  z=Lice.boot[i,,j,2], pmat=pmat1)
	Lmax<-trans3d(x=dist, y=day,  z=Lice.boot[i,,j,3], pmat=pmat1)
	
	for(k in 1:length(day)){
		ind<-c(findInterval(dist[k],seq(-44,68,d)), findInterval(day[k],seq(100,144,d)))
		
		# text(Lmean$x[k], Lmean$y[k], k, xpd=NA)
			
		if(Lice.boot[i,k,j,1]>=Z[ind[1], ind[2]]){ #if obs>pred
			points(Lmean$x[k], Lmean$y[k], pch=19, xpd=NA)
			segments(Lmean$x[k], Lmean$y[k], Lmax$x[k], Lmax$y[k])
			
			if(Lice.boot[i,k,j,2]>Z[ind[1], ind[2]]){
				segments(Lmean$x[k], Lmean$y[k], Lmin$x[k], Lmin$y[k])
			}else{
				modMin<-trans3d(x=dist[k], y=day[k], z=Z[ind[1], ind[2]], pmat=pmat1)
				segments(Lmean$x[k], Lmean$y[k], modMin$x, modMin$y)
				segments(modMin$x, modMin$y, Lmin$x[k], Lmin$y[k], col="#00000040" )
				}
			}
		if(Lice.boot[i,k,j,1]<Z[ind[1], ind[2]]){ #if obs>pred	
			points(Lmean$x[k], Lmean$y[k], pch=19, col="#00000040")
			segments(Lmean$x[k], Lmean$y[k], Lmin$x[k], Lmin$y[k], col="#00000040")
			
			if(Lice.boot[i,k,j,3]<Z[ind[1], ind[2]]){
				segments(Lmean$x[k], Lmean$y[k], Lmax$x[k], Lmax$y[k], col="#00000040")
			}else{
				modMin<-trans3d(x=dist[k], y=day[k], z=Z[ind[1], ind[2]], pmat=pmat1)
				segments(Lmean$x[k], Lmean$y[k], modMin$x, modMin$y, col="#00000040")
				segments(modMin$x, modMin$y, Lmax$x[k], Lmax$y[k] )
				}
			}
			
		}
	mtext(side=3, adj=0, c("a) Copepodid", "b) Chalimus", "c) Motile")[i])	
	
} #end stage i

#par(new=TRUE)
#persp(z=matrix(sim.out[[i]][,j], nrow=n.dist, ncol=n.day, byrow=FALSE), x=seq(-44,68,d), y=seq(100,144,d), xlab="", ylab="", zlab="", theta=-150, phi=25, col="#FFFFFF80", border=NA, bty="n")

# source('~/Desktop/filledContour.R')
# greys<-colorRampPalette(c("white", 1))
# filledContour(z=matrix(sim.out[[stage]][,species], nrow=n.dist, ncol=n.day, byrow=FALSE), x=seq(-44,68,d), y=seq(100,144,d), xlab="Distance", ylab="Time", zlab=c("C(x,t)", "H(x,t)", "M(x,t)")[stage], main=c("Pink", "Chum")[species], color.palette=greys, nlevels=20)
# contour(z=matrix(sim.out[[stage]][,species], nrow=n.dist, ncol=n.day, byrow=FALSE), x=seq(-44,68,d), y=seq(100,144,d), add=TRUE)

#######################################################
# Animate model results
#######################################################
library(animation)

dist2<-rep(c(-44:68), length(unique(day)))
day2<-rep(sort(unique(day)), each=113)
sim.ani<-simulate.lice(p=S[[1]][,1], dist=dist2, day=day2)

ani.options(ani.width=4.5, ani.height=6, ani.dev="pdf", ani.type="pdf")
setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model/Results/Animation")
saveLatex(expr={
		# Lice.boot[stage, site.id, species, metric(mean, lci, uci)]
for(i in 1:length(unique(day))){
	par(mfrow=c(3,1), mar=c(2,4,2,1), oma=c(3,1,2,0))
	siteId<-which(day==sort(unique(day))[i])
	for(s in 1:3){
		plotCI(dist[siteId], Lice.boot[s,siteId,1,1], li=Lice.boot[s,siteId,1,2], ui=Lice.boot[s,siteId,1,3], xlim=range(dist), ylim=c(0, max(Lice.boot[s,,,])), col="#FF000080", gap=0.3, xlab="", ylab=c("C(x,t)", "H(x,t)", "M(x,t)")[s], las=1, bty="l", cex.lab=1.5)
		plotCI(dist[siteId], Lice.boot[s,siteId,2,1], li=Lice.boot[s,siteId,2,2], ui=Lice.boot[s,siteId,2,3], add=TRUE, col="#0000FF80", gap=0.3)
		lines(dist2[which(day2==sort(unique(day))[i])], sim.ani[[s]][which(day2==sort(unique(day))[i]),1], col="#FF000080", lwd=2)
		lines(dist2[which(day2==sort(unique(day))[i])], sim.ani[[s]][which(day2==sort(unique(day))[i]),2], col="#0000FF80", lwd=2)
		mtext(side=3, adj=0, paste("  ", letters[s], ")", sep=""))
		if(s==1) legend("topleft", lwd=2, c("Pink", "Chum"), col=c("#FF000080", "#0000FF80"), bty="n")
		} #end s stage
	mtext(side=1, "Distance (km)", line=3)
	mtext(side=3, outer=TRUE, strftime(as.Date(sort(unique(day))[i], origin="2006-01-01"), format="%B %d, %Y"))
	} #END I
},pdflatex = '/usr/texbin/pdflatex', img.name="GreensFit", latex.filename="GreensFit_ani.tex", interval=0.3, caption="Observations (points; $\\pm$ 95\\% CI) and model predictions (lines) for the number of (a) copepodid, b) chalimus and c) motile lice per wild juvenile pink salmon (red line) and chum salmon (blue line).")
	
