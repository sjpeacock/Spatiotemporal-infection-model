rm(list=ls())
# setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model")
setwd("~/Greens/Spatiotemporal-infection-model")

library(dclone)
library(parallel)

source("0 Original model/model_original2.R")

########################################################
# 1) Read in data and define parameters
########################################################

#--------------------------------------------------------------------------------------
# Read in index file with space and time for sites
index = read.csv("3 Fitting/index.csv", header=TRUE);
day<-index$day; dist<-index$dist

#--------------------------------------------------------------------------------------
# Read in Lice data
data <-read.csv("3 Fitting/lice.csv", header=TRUE, na.string="NA")
n.obs <- length(data$C)
n.sites <- max(data$site.id)
n.reps <- length(unique(index$rep))

# data$rep<-c(
	# rep(1, length(which(data$site<=max(index$site.id[index$rep==1])))),
	# rep(2, length(which(data$site>max(index$site.id[index$rep==1])&data$site<=max(index$site.id[index$rep==2])))), 
	# rep(3, length(which(data$site>max(index$site.id[index$rep==2])&data$site<=max(index$site.id[index$rep==3]))))) 
	
#########################################################
# 2) Run Bayesian model fitting
#########################################################
#Priors & RNG
set.seed(314159)
mu<-c(1,1,0,0,1,1,0,0,1); sig<-c(1,1,2,1.7,1.7,1.7,1.7,1.7,1.1)
par.names<-c("k.log", "phi.log", "sc.logit", "sh.logit", "Lc.log", "D.log", "Lh_R.log", "Lm_R.log", "r.log")
init.p<-list(); length(init.p)<-3
# Start with initial values overdispersed w.r.t. prior distribtions
for(i in 1:3){
	init.p[[i]]<-list();
	init.p[[i]]$k.log<-rnorm(2, mu[1], sig[1]+0.1)
	init.p[[i]]$phi.log<-matrix(rnorm(9, mu[2], sig[2]+0.1), nrow=3, ncol=3)
	init.p[[i]]$sc.logit<-rnorm(2, mu[4], sig[4]+0.1)
	init.p[[i]]$sh.logit<-rnorm(2, mu[5], sig[5]+0.1)
	init.p[[i]]$D.log<-rnorm(1, mu[6], sig[6]+0.1)
	init.p[[i]]$Lh_R.log<-rnorm(1, mu[7], sig[7]+0.1)
	init.p[[i]]$Lm_R.log<-rnorm(1, mu[8], sig[8]+0.1)
	init.p[[i]]$r.log<-rnorm(1, mu[9], sig[9]+0.1)
	
	#init.p[[i]]$.RBG.name<-c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper")[i]

	#init.p[[i]]$.RNG.seed<-list(c(13129,  2996,  9875), c(1387288580, 1116634134), c(790609190, 186962253))[[i]]
		
	}
# par(mfrow=c(2,4)); x.max<-c(50,50,50,1,1,1,1,40)
# for(i in c(1:3,8)) hist(exp(rnorm(10^4, mu[i], sig[i])), main=paste(i, par.names[i]), xlab="", breaks=50, xlim=c(0,x.max[i]))	
# for(i in 4:7)hist(inv.logit(rnorm(10^4, mu[i], sig[i])), main=paste(i, par.names[i]), xlab="", breaks=50)	
# Use initial values that worked for the shorter run!

dat<-list(y=c(-3.7, 4, 53), 
	gamma=1.56, 
	u_n=4/5, 
	u_c=1/5, 
	x=dist, 
	n.sites=n.sites,
	n.reps=n.reps,
	C=data$C, 
	H=data$H,  
	M=data$M,
	species=data$spp, 
	site=data$site, 
	rep=index$rep,
	mu=mu, tau=sig^(-2)) 

cl<-makeCluster(6)
t.start<-proc.time()[3]

#mod<-jags.parfit(cl, data=dat, params=c("k.log", "phi.log", "sc.logit", "sh.logit", "Lc.log", "D.log", "Lh_R.log", "Lm_R.log", "r.log"),  model=model, inits=init.p, n.adapt=1000, n.update=1000, n.iter=1000, n.chains=3)

mod<-dc.parfit(cl, data=dat, params=c("k.log", "phi.log", "sc.logit", "sh.logit", "Lc.log", "D.log", "Lh_R.log", "Lm_R.log", "r.log"),  model=model, inits=init.p, n.adapt=5000, n.update=10000, n.iter=5000, n.chains=3, unchanged=c("n.sites", "n.reps", "gamma", "u_n", "u_c", "mu", "tau", "y"), n.clones=c(1,2,3,4,5,6,8,10))

cat(paste("Process time (minutes) = ", round((proc.time()[3]-t.start)/60, 4)))

stopCluster(cl)

save.image("OriginalModel_DClone_20161219.RData")


#########################################
