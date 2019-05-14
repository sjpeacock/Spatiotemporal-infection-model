### generate the distributions of copepods from the salmon farms in 2006
#-------------------------------------------------------------------------------------------------
library(snowfall)
library(lattice)

source("2 Larval-distribution/functions.R")
source("grid.R")
source("params.R")
rm(p.lice)
#------------------------------------------------------------------------------------------------
#Number of cores to use in the cope caluculation
cpus<-4
 
###############################################################################
## Create copepodid distribution in space and time for different farms
###############################################################################
# ***** time = 0 is January 1, 2006 *******
# ***** space = 0 is Shewell Island (used to be SE Viscount but I changed this in June 2011) *******

# This takes ~10 minutes per farm using 4 cores

#------------------------------------------------------------------------------------------------
# Sargeaunts Passage farm

Xsource = -3.71		# location of farm 
Tsource = as.numeric(as.Date("2006-03-07")-as.Date("2006-01-01"))  # treatment date: March 7, 2006 
zeroindex = which(T==Tsource)
r1 = p.farm[[1]][1,1]	# rate of exponential growth
r2 = p.farm[[1]][2,1]	# rate of exponential decay
f0 = p.farm[[1]][3,1]	# source strength at time 0

# !This takes ~ 10 minutes !
t.start<-proc.time()
source("2 Larval-distribution/sim_larvae_1farm.R")
cat("Process time for whole SP copepodid distribution (minutes) = ", (proc.time()-t.start)[3]/60)

SP = est;
SPmaxarea = sum(est[,zeroindex])*dx;

#wireframe(SP, row.values=x, column.values=T, col="#FF000020", scales=list(arrows=FALSE))

#-------------------------------------------------------------------------------------------------

# Humphrey Rock farm
Xsource = 4		# location of farm
Tsource = as.numeric(as.Date("2006-02-13")-as.Date("2006-01-01")) # treatment date in days (Feb. 13, 2006 (zero is Feb 17, 2006))
zeroindex = which(T==Tsource)
r1 = p.farm[[2]][1,1]	# rate of exponential growth
r2 = p.farm[[2]][2,1];	# rate of exponential decay

f0 = p.farm[[2]][3,1]; #source strenght at time 0

# !This takes ~ 10 minutes !
t.start<-proc.time()
source("2 Larval-distribution/sim_larvae_1farm.R")
cat("Process time for whole HR copepodid distribution (minutes) = ", (proc.time()-t.start)[3]/60)

HR = est;
HRmaxarea = sum(est[,zeroindex])*dx;

#------------------------------------------------------------------------------------------------------------

# Burdwoods farm - no data
Xsource = 53;	# location of farm
Tsource = as.numeric(as.Date("2006-01-13")-as.Date("2006-01-01"));	# treatment date in days (January 13, 2006)
zeroindex = which(T==Tsource)
r1 = p.farm[[3]][1,1];	# rate of exponential growth(estimated from Cohen Commission data)
r2 = p.farm[[3]][2,1];	# rate of exponential decay

f0 = p.farm[[3]][3,1]; #source strenght at time 0

# !This takes ~ 10 minutes !
t.start<-proc.time()
source("2 Larval-distribution/sim_larvae_1farm.R")
cat("Process time for whole BW copepodid distribution (minutes) = ", (proc.time()-t.start)[3]/60)

BW = est;
BWmaxarea = sum(est[,zeroindex])*dx;

#------------------------------------------------------------------------------------------------------------
# Write to .csv file
# Scale output according the whichever farm had more lice
scale.par<-max(c(SPmaxarea, HRmaxarea, BWmaxarea))

copedist<-SP/scale.par+HR/scale.par+BW/scale.par
write.csv(copedist, paste("2 Larval-distribution/allFarms06_", Sys.Date(), "=.csv", sep=""), row.names=FALSE);


#write.csv(SP/scale.par, "Output/SP06.csv", row.names=FALSE, col.names=FALSE);
#write.csv(HR/scale.par, "Output/HR06.csv", row.names=FALSE, col.names=FALSE);
#write.csv(BW/scale.par, "Output/BW06.csv", row.names=FALSE, col.names=FALSE);
#
##------------------------------------------------------------------------------------------------------------
## Plot distributions
#
#wireframe(copedist, row.values=x, column.values=T, col="#00000010", scales=list(arrows=FALSE))
#
#
#SP.subset<-SP[seq(1, 4096, 40), seq(1, 1600, 20)]/HRmaxarea
#wireframe(SP.subset, row.values=seq(1, 4096, 40), col.values=seq(1, 1600, 20), xlab="x", ylab="t", zlab="larvae", main="Sargeaunts")
#
#HR.subset<-HR[seq(1, 4096, 40), seq(1, 1600, 20)]/HRmaxarea
#wireframe(HR.subset, row.values=seq(1, 4096, 40), col.values=seq(1, 1600, 20), xlab="x", ylab="t", zlab="larvae",main="Humphrey")
#
#allF<-SP.subset+HR.subset
#wireframe(copedist[seq(1, 4096, 40), seq(1, 1600, 20)], row.values=seq(1, 4096, 40), col.values=seq(1, 1600, 20), xlab="x", ylab="t", zlab="larvae", main="Both farms")