setwd("~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model/Farm model")

lice_data<-read.csv("farm_lice.csv")
lice_data$date<-as.Date(paste(lice_data$year, lice_data$month, lice_data$day, sep="-"))

treatment.date<-c(as.Date("2006-03-07"), as.Date("2006-02-13"), as.Date("2006-01-13"))

library(dclone)


model<-function(){
	# Predicted values
	for(i in 1:length(sampling.date)){
   		f[i]<-f0*(ind[i]*exp(r[1]*(sampling.date[i]-treatment.date))+(1-ind[i])*exp(r[2]*(sampling.date[i]-treatment.date)))
   		}  	
   	# Likelihood calculation
   	f.prob <- k/(k+n.fish*f)
   	for(i in 1:length(y)){
   		y[i] ~ dnbinom(f.prob[sampling.date.ind[i]], k)
   	}
   	# Priors on parameters
	for(i in 1:2){r[i] ~ dnorm(0, 0.1)} # growth rate
	f0 ~ dlnorm(0, 0.1) # maximum louse abundance at treatment
	k ~ dlnorm(0, 0.1) #overdispersion parameter in negative binomial
}

#j = farm number 1-3 where 1=SP, 2=HR, 3=BG
allFits<-list(); length(allFits)<-3
for(j in 1:2){
	
	farm.ind<-which(lice_data$farm==c("SP", "HR", "BG")[j])
	dat<-list(
		sampling.date=unique(as.numeric(lice_data$date[farm.ind])), 
		treatment.date=as.numeric(treatment.date[j]), 
		y=lice_data$lice[farm.ind], 
		n.fish=unique(lice_data$num_fish[farm.ind]))
	
	dat$ind<-as.numeric(dat$sampling.date<dat$treatment.date)
	
	dat$sampling.date.ind<-match(as.numeric(lice_data$date[farm.ind]), dat$sampling.date)
	
	cl<-makeCluster(4, type="SOCK")
	t.start<-proc.time()

	dc.out<-dc.parfit(cl, data=dat, params=c("r", "f0", "k"), model=model, n.clones=c(1:10), n.chains=3, n.adapt=1000, n.update=1000, n.iter=5000, unchanged=c("sampling.date", "treatment.date", "n.fish", "ind"))
	
	all.time<-(proc.time()-t.start)[3]/60
	cat("Process time (minutes) = ", all.time)
	stopCluster(cl)

	allFits[[j]]<-dc.out
	
	}

##############################################
# Different approach for BG
##############################################

model.norm<-function(){
	# Predicted values
	for(i in 1:length(sampling.date)){
   		f[i]<-f0*(ind[i]*exp(r[1]*(sampling.date[i]-treatment.date))+(1-ind[i])*exp(r[2]*(sampling.date[i]-treatment.date)))
   		}  	
   	# Likelihood calculation
   	for(i in 1:length(y)){
   		y[i] ~ dnorm(f[sampling.date.ind[i]], pow(k, -2))
   	}
   	# Priors on parameters
	for(i in 1:2){r[i] ~ dnorm(0, 0.1)} # growth rate
	f0 ~ dlnorm(0, 0.1) # maximum louse abundance at treatment
	k ~ dlnorm(-1, 0.1) #standard deviation in normal distribution
}
j<-3
farm.ind<-which(lice_data$farm=="BG")
dat<-list(
		sampling.date=unique(as.numeric(lice_data$date[farm.ind])), 
		treatment.date=as.numeric(treatment.date[j]), 
		y=lice_data$lice[farm.ind]/lice_data$num_fish[farm.ind])
dat$ind<-as.numeric(dat$sampling.date<dat$treatment.date)
dat$sampling.date.ind<-match(as.numeric(lice_data$date[farm.ind]), dat$sampling.date)

cl<-makeCluster(4, type="SOCK")
t.start<-proc.time()

dc.out<-dc.parfit(cl, data=dat, params=c("r", "f0", "k"), model=model.norm, n.clones=c(1:10), n.chains=3, n.adapt=1000, n.update=1000, n.iter=5000, unchanged=c("sampling.date", "treatment.date", "n.fish", "ind"))

all.time<-(proc.time()-t.start)[3]/60
cat("Process time (minutes) = ", all.time)
stopCluster(cl)

allFits[[3]]<-dc.out
	
##############################################
# print parameters
##############################################
p.farm<-list(); length(p.farm)<-3
pOrd<-c(3,4,1,2)

for(i in 1:3){
	S<-summary(allFits[[i]])[[1]]
	p.farm[[i]]<-cbind(S[pOrd,1], S[pOrd,1]-1.96*S[pOrd,3], S[pOrd,1]+1.96*S[pOrd,3])
	rownames(p.farm[[i]])<-rownames(S)[pOrd]
	colnames(p.farm[[i]])<-c("Est", "2.5%", "97.5%")
	p.farm[[i]]<-round(p.farm[[i]], 3)
}

dump("p.farm")
