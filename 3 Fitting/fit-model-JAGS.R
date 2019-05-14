library(here)
library(dclone)
library(parallel)

source("3 Fitting/model.R")

########################################################
# 1) Read in data and define parameters
########################################################

#--------------------------------------------------------------------------------------
# Read in index file with space and time for sites
index <- read.csv("3 Fitting/index.csv", header=TRUE)
day <- index$day
dist <- index$dist

#--------------------------------------------------------------------------------------
# Read in Lice data
data <-read.csv("3 Fitting/lice.csv", header=TRUE, na.string="NA")
n.obs <- length(data$C)
n.sites <- max(data$site.id)

#--------------------------------------------------------------------------------------
# Define fixed parameters
source("grid.R")

#--------------------------------------------------------------------------------------
# Read in calculated larval distributions
# "_trunc" version just has fewer digits (= smaller file), but same distributions
farmL <- read.csv("2 Larval-distribution/larvalDistribution2006_trunc.csv", header=FALSE) 
# farmL <- read.csv("2 Larval-distribution/larvalDistribution2006_2019-01-08.csv", header=FALSE) 
farmL <- as.matrix(farmL)
dimnames(farmL)[[2]] <- NULL
# Add a row of zeroes for summing purposes.
farmL <- rbind(rep(0, length(T)), farmL)

#Fixed parameters from Stien et al. (2005)
tauC <- 0.525^-2 # Development time of copepodids
tauH <- 0.250^-2 # Development time of chalimus
tauM <- 0.187^-2 # Development time of motiles

#Indices for summing over time (do not depend on free parameters)
C_index <- c(0:(round(tauC / dt) - 1))
H_index <- c(round(tauC / dt):(round(tauH / dt) - 1))
M_index <- c(round(tauH / dt):(round(tauM / dt) - 1))

#########################################################
# 2) Run Bayesian model fitting
#########################################################
#Priors & RNG
set.seed(314159)
mu <- c(1,1,0,0,0,0,0,1) # mean in prior distributions
sig<-c(1,1,2,1.7,1.7,1.7,1.7,1.1) #sd in prior distributions
par.names<-c("k.log", "phi.log", "V.log", "Psc.log", "Psh.log", "Csc.log", "Csh.log", "r.log")
init.p<-list(); length(init.p)<-3

# Start with initial values overdispersed w.r.t. prior distribtions
for(i in 1:3){
	init.p[[i]]<-list();
	init.p[[i]]$k.log<-rnorm(1, mu[1], sig[1]+0.1)
	init.p[[i]]$phi.log<-rnorm(1, mu[2], sig[2]+0.1)
	init.p[[i]]$V.log<-rnorm(1, mu[3], sig[3]+0.1)
	init.p[[i]]$Psc.log<-rnorm(1, mu[4], sig[4]+0.1)
	init.p[[i]]$Psh.log<-rnorm(1, mu[5], sig[5]+0.1)
	init.p[[i]]$Csc.log<-rnorm(1, mu[6], sig[6]+0.1)
	init.p[[i]]$Csh.log<-rnorm(1, mu[7], sig[7]+0.1)
	init.p[[i]]$r.log<-rnorm(1, mu[8], sig[8]+0.1)
	
	# Uncomment to set random number sequence 
	# Don't want this for data cloning...
	#init.p[[i]]$.RBG.name<-c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper")[i]
	#init.p[[i]]$.RNG.seed<-list(c(13129,  2996,  9875), c(1387288580, 1116634134), c(790609190, 186962253))[[i]]
		
	}

# Plot priors
# par(mfrow=c(2,4)); x.max<-c(50,50,50,1,1,1,1,40)
# for(i in c(1:3,8)) hist(exp(rnorm(10^4, mu[i], sig[i])), main=paste(i, par.names[i]), xlab="", breaks=50, xlim=c(0,x.max[i]))	
# for(i in 4:7)hist(inv.logit(rnorm(10^4, mu[i], sig[i])), main=paste(i, par.names[i]), xlab="", breaks=50)	

#--------------------------------------------------------------------------------------
# Data list to be passed to JAGS:
dat<-list(dt=dt, dx=dx, 
	L=farmL, 
	site.id=data$site.id, spp=data$spp,
	n.sites=n.sites,
	C=data$C, H=data$H, M=data$M,
	dist=dist, day=day,
	x0=which(x==0), t0=0,
	C_index=c(0:(round(tauC/dt)-1)), 
	H_index=c(round(tauC/dt):(round(tauH/dt)-1)), 
	M_index=c(round(tauH/dt):(round(tauM/dt)-1)),
	tauC=tauC, tauH=tauH, tauM=tauM,
	mu=mu, tau=sig^(-2)) 

#--------------------------------------------------------------------------------------
# Model fitting
cl<-makeCluster(4)
t.start<-proc.time()[3]

mod<-dc.parfit(cl, data=dat, params=c("phi.log", "V.log", "k.log", "Psc.log", "Psh.log", "Csc.log", "Csh.log", "r.log"),  model=model, inits=init.p, n.adapt=5000, n.update=20000, n.iter=2000, n.chains=3, unchanged=c("dt", "dx", "L", "x0", "t0", "tauC", "tauH", "tauM", "C_index", "H_index", "M_index", "dist", "day", "n.sites", "mu", "tau"), n.clones=c(1:6, 8, 10))

ptime <- round((proc.time()[3]-t.start)/60
cat(paste("Process time (minutes) = ", ptime, 4)))

stopCluster(cl)

# Takes approximately 7 days on 4 processors, with K = 1-6, 8, 10

# Save MCMC output for "Looking-at-JAGS-output.R"
save.image(paste("LiceFits2Data_", Sys.Date(), ".RData", sep=""))


