######################################################################################
## Figures for Green's function model paper
######################################################################################

#-------------------------------------------------------------------------------------
## Figure 1: Map
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
## Figure 1: Model schematic
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
## Figure 3: Migration path
#-------------------------------------------------------------------------------------
load("~/Google Drive/Greens/Code/Simulations/FarmLocSimulations_2012Aug22.RData")

#days migrating = 80 km/v
migration.time<-floor(80/p.lice[1,5])

t.migrating<-seq(73, (73+migration.time), dt)
as.Date(t.migrating, origin=c("2006-01-01"))
x.migrating<-round(-20+(t.migrating-73)*p.lice[1,5], 1)

ind.migrating.attached<-cbind(findInterval(x.migrating, x.sim0), findInterval(t.migrating, t.sim0))
ind.migrating<-cbind(findInterval(x.migrating, x), findInterval(t.migrating, T))

x.migrating<-x.migrating+20
x.sim0<-x.sim0+20
x<-x+20
#----------
# Function
plot3D<-function(Z, zmax){
	z   <- Z$z
	nrz <- nrow(z)
	ncz <- ncol(z)
	zlim=c(0,zmax)
	Ncol<-10

	couleurs  <- paste(tail(grey(seq(0, trunc(1.2 * Ncol), 1)/(1.2*Ncol)),Ncol), rep("80", 1.2*Ncol), sep="") #transparent greyscale
	fcol      <- couleurs[as.matrix(trunc(z/zlim[2]*(Ncol-1))+1)]
	dim(fcol) <- c(nrz,ncz)
	fcol      <- fcol[-nrz,-ncz]
	lcol<-1
	
	pmat1<-persp(Z,col=fcol,zlim=zlim,theta=135,phi=30,zlab="", xlab="", ylab="", ticktype="detailed", border=lcol, cex.axis=0.8)
	text(0.27, -0.43, srt=46, "Distance (km)")
	text(-0.3, -0.45, srt=-46, "Time (days)")
	return(pmat1)
	}

#----------
## Plots of migration path

quartz(width=4, height=7)

par(mfrow=c(2,1), mar=c(2, 3, 1, 0.5), oma=rep(0,4), mai=c(0.4, 0.2, 0.2, 0.1))
thin<-15
f=1;m=M; s=1
Z1<-FreeCope[[f,m]]/scale.par
Z<-list(z=100*Z1[seq(1, length(x), thin), seq(1, length(T), thin)], x=x[seq(1, length(x), thin)], y=T[seq(1, length(T), thin)])
pmat1<-plot3D(Z, zmax=1)
text(-0.52, -0.03, srt=102, "Infectious copepodites")
lines(trans3d(x=x.migrating, y=t.migrating, z=100*Z1[ind.migrating], pmat=pmat1), col=2)
start.<-ind.migrating[1,]; dim(start.)<-c(1,2)
points(trans3d(x=x.migrating[1], y=t.migrating[1], z=Z1[start.], pmat=pmat1), col=2, pch=21, bg="white", cex=0.8)
end<-ind.migrating[length(x.migrating),]; dim(end)<-c(1,2)
points(trans3d(x=x.migrating[length(x.migrating)], y=t.migrating[length(x.migrating)], z=Z1[end], pmat=pmat1), col=2, pch=19, cex=0.8)
mtext(side=3, adj=0, "(a)")


thin<-10
Z1<-matrix(Attached[[f,m]][[3]][,s], nrow=length(x.sim0), ncol=length(t.sim0), byrow=TRUE)
Z<-list(z=Z1[seq(1, length(x.sim0), thin), seq(1, length(t.sim0), thin)], x=x.sim0[seq(1, length(x.sim0), thin)], y=t.sim0[seq(1, length(t.sim0), thin)])
pmat1<-plot3D(Z, zmax=0.5)
text(-0.52, -0.03, srt=102, "Attached motiles")
lines(trans3d(x=x.migrating, y=t.migrating, z=Z1[ind.migrating.attached], pmat=pmat1), col=2)
start.<-ind.migrating.attached[1,]; dim(start.)<-c(1,2)
points(trans3d(x=x.migrating[1], y=t.migrating[1], z=Z1[start.], pmat=pmat1), col=2, pch=21, bg="white", cex=0.8)
end<-ind.migrating.attached[length(x.migrating),]; dim(end)<-c(1,2)
points(trans3d(x=x.migrating[length(x.migrating)], y=t.migrating[length(x.migrating)], z=Z1[end], pmat=pmat1), col=2, pch=19, cex=0.8)
mtext(side=3, adj=0, "(b)")
#----------


#-------------------------------------------------------------------------------------
## Figure 4: Farm predictions
#-------------------------------------------------------------------------------------
require(gplots)
load("~/Google Drive/Greens/FINALLY/Figures/farmModels_13June2012.RData")
quartz(width=4, height=7)
par(mfrow=c(3,1), mar=c(4,3,1,1), oma=c(2,2,1,1))

plot(as.Date(T.all, origin="1970-01-01"), fset(p.fitted[[1]][,1], T.all, T.slice.all[1]), "l", xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2)
polygon(x=c(T.all, rev(T.all)), y=c(CI.pred[[1]][,1], rev(CI.pred[[1]][,2])), border=NA, col="#00000030")
plotCI(as.Date(unique(Z$Date[Z$Farm=="SP"]), origin="1970-01-01"), Lice.summary[[1]][,1], ui=Lice.summary[[1]][,3], li=Lice.summary[[1]][,2], add=TRUE)
abline(v=T.slice.all[1], lty=2)
mtext(side=3, adj=0, "(a) Sargeaunt's Passage", line=0.5)

plot(as.Date(T.all, origin="1970-01-01"), fset(p.fitted[[2]][,1], T.all, T.slice.all[2]), "l", xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2)
polygon(x=c(T.all, rev(T.all)), y=c(CI.pred[[2]][,1], rev(CI.pred[[2]][,2])), border=NA, col="#00000030")
plotCI(as.Date(unique(Z$Date[Z$Farm=="HR"]), origin="1970-01-01"), Lice.summary[[2]][,1], ui=Lice.summary[[2]][,3], li=Lice.summary[[2]][,2], add=TRUE)
abline(v=T.slice.all[2], lty=2)
mtext(side=3, adj=0, "(b) Humphrey Rock", line=0.5)

plot(as.Date(T.all, origin="1970-01-01"), fset(p.fitted[[3]][,1], T.all, T.slice.all[3]), "l", xlab="", ylab="", bty="n", ylim=c(0,8), las=1, cex.axis=1.2)
polygon(x=c(T.all, rev(T.all)), y=c(CI.pred[[3]][,1], rev(CI.pred[[3]][,2])), border=NA, col="#00000030")
points(as.Date(Z$Date[Z$Farm=="BG"], origin="1970-01-01"), Z$Lice[Z$Farm=="BG"]/20)
abline(v=T.slice.all[3], lty=2)
mtext(side=3, adj=0, "(c) Burdwood Group", line=0.5)

mtext(side=2, outer=TRUE, "Average lice per farm salmon")
mtext(side=1, outer=TRUE, "Date (2005/2006)")


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

