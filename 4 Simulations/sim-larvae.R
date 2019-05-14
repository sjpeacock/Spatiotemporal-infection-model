### generate the distributions of copepods from the salmon farms in 2006
#-------------------------------------------------------------------------------------------------
setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model")
source("4 Simulations/sim-functions.R")
source("4 Simulations/sim-grid.R")
library(snowfall)
library(lattice)

source("params.R")
#------------------------------------------------------------------------------------------------
#Number of cores to use in the cope caluculation
cpus<-4
 
 farmL.all<-list(); length(farmL.all)<-4
#------------------------------------------------------------------------------------------------
for(s in 1:4){# for each scenario 
	#---------------------
	# Sargeaunts Passage farm
	#---------------------
	Xsource = -3.71		# location of farm 
	t_treat = round(tT[[s]][1])  # treatment date 
	zeroindex = which(T==t_treat)
	r1 = p.farm[[1]][1,1]	# rate of exponential growth
	r2 = p.farm[[1]][2,1]	# rate of exponential decay
	f_treat = fT[[s]][1]			# source strength at time 0
	
	# !This takes ~ 10 minutes !
	t.start<-proc.time()
	source("4 Simulations/sim-larvae_1farm.R")
	cat("Process time for whole SP copepodid distribution (minutes) = ", (proc.time()-t.start)[3]/60)
	
	SP = est;
	SPmaxarea = sum(est[,zeroindex])*dx;
	
	#---------------------
	# Humphrey Rock farm
	#---------------------
	Xsource = 4		# location of farm
	t_treat = round(tT[[s]][2])  # treatment date 
	zeroindex = which(T==t_treat)
	r1 = p.farm[[2]][1,1]	# rate of exponential growth
	r2 = p.farm[[2]][2,1]	# rate of exponential decay
	f_treat = fT[[s]][2]			# source strength at time 0
	
	# !This takes ~ 10 minutes !
	t.start<-proc.time()
	source("4 Simulations/sim-larvae_1farm.R")
	cat("Process time for whole HR copepodid distribution (minutes) = ", (proc.time()-t.start)[3]/60)
	
	HR = est;
	HRmaxarea = sum(est[,zeroindex])*dx;
	
	#---------------------
	# Burdwoods farm 
	#---------------------
	Xsource = 53;	# location of farm
	t_treat = round(tT[[s]][3])  # treatment date 
	zeroindex = which(T==t_treat)
	r1 = p.farm[[3]][1,1]	# rate of exponential growth
	r2 = p.farm[[3]][2,1]	# rate of exponential decay
	f_treat = fT[[s]][3]			# source strength at time 0
	
	# !This takes ~ 10 minutes !
	t.start<-proc.time()
	source("4 Simulations/sim-larvae_1farm.R")
	cat("Process time for whole BW copepodid distribution (minutes) = ", (proc.time()-t.start)[3]/60)
	
	BW = est;
	BWmaxarea = sum(est[,zeroindex])*dx;

	#---------------------------------------------------------------------------------------------------
	# Write to .csv file
	# Scale output according the whichever farm had more lice
	scale.par<-max(c(SPmaxarea, HRmaxarea, BWmaxarea))

	copedist<-SP/scale.par+HR/scale.par+BW/scale.par
	
	farmL.all[[s]]<-copedist
	} #end scenario