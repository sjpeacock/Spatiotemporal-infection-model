# This copy attemps to fit all three reps, using most of the same parameters but allowing the farm source to change between replicates
sim.model<-function(parms, dat){
	#Unlist data
	y<-dat$y
	gamma<-dat$gamma
	u_n<-dat$u_n
	u_c<-dat$u_c
	x<-dat$x
	n.sites<-dat$n.sites
	n.reps<-dat$n.reps
	# species<-dat$species
	# site<-dat$site
	# rep<-dat$rep
	
	#-----------------------------------------------------------------------------------------------	
	# Parameters
	# Species specific
	k<-c(exp(parms['k.log[1]']), exp(parms['k.log[2]']))
	sc<- c(exp(parms['sc.logit[1]'])/(1+exp(parms['sc.logit[1]'])), exp(parms['sc.logit[2]'])/(1+exp(parms['sc.logit[2]'])))
	sh <- c(exp(parms['sh.logit[1]'])/(1+exp(parms['sh.logit[1]'])), exp(parms['sh.logit[2]'])/(1+exp(parms['sh.logit[2]'])))
	Lc<-c(exp(parms['Lc.log[1]']), exp(parms['Lc.log[2]']))
		
	# Global parameters
	phi<-matrix(nrow=n.reps, ncol=length(y))
	for(i in 1:n.reps){
		for(j in 1:length(y)){ #for each farm
			phi[i,j] <- exp(parms[which(names(parms)==paste("phi.log[", i, ",",j, "]", sep=""))])
		}}
	
	D<-exp(parms['D.log'])
	Lh_R<-exp(parms['Lh_R.log'])
	Lm_R<-exp(parms['Lm_R.log'])
	r<-exp(parms['r.log'])

	# Define parameters for larval distribution
	a1<-(gamma+(gamma^2+4*u_n*D)^0.5)/(2*D); a2<-abs((gamma-(gamma^2+4*u_n*D)^0.5)/(2*D));
	b1<-(gamma+(gamma^2+4*u_c*D)^0.5)/(2*D); b2<-abs((gamma-(gamma^2+4*u_c*D)^0.5)/(2*D));
	
	# Larval constant (integral from -infty to infty)
	c1<-1/(a1*(a1+b2))+1/(b1*(a1-b1))-1/(a1*(a1-b1))+1/(b1*(a2+b1))+1/(b2*(a1+b2))+1/(b2*(a2-b2))-1/(a2*(a2-b2))+1/(a2*(a2+b1))
	
	#----------------------------------------------------------------------------------------
## Farm Footprints
# Create indices of which solution to use (left hand, right hand, or middle)
# The solution is in three parts, and we calculate all three and then the index tells us which of the parts we keep.
 index<-matrix();length(index)<-3*n.sites*3*length(y)*2;dim(index)<-c(3,n.sites, 3, length(y),2)
for(j in 1:length(y)){ # for each farm
	for(i in 1:n.sites){ # for each site		
		for(s in 1:2){ # for each species (1=pink, 2=chum)
			# Copepidite
			index[1,i,1,j,s]<-(x[i]<=y[j])
			index[1,i,2,j,s]<-((x[i]-Lc[s]) <= y[j]) && (x[i]>y[j])
			index[1,i,3,j,s]<-((x[i]-Lc[s]) > y[j])
			# Chalimus
			index[2,i,1,j,s]<-((x[i]-Lc[s]) <= y[j])
			index[2,i,2,j,s]<-((x[i]-Lc[s]*(1+Lh_R)) <= y[j] && (x[i]-Lc[s]) > y[j])
			index[2,i,3,j,s]<-((x[i]-Lc[s]*(1+Lh_R)) > y[j])
			# Motile
			index[3,i,1,j,s]<-((x[i]-Lc[s]*(1+Lh_R)) <= y[j])
			index[3,i,2,j,s]<-((x[i]-(Lc[s]*(1+Lh_R+Lm_R))) <= y[j] && (x[i]-Lc[s]*(1+Lh_R)) > y[j])
			index[3,i,3,j,s]<-((x[i]-(Lc[s]*(1+Lh_R+Lm_R))) > y[j])
			} # end s species	
			}# end i sites
		}# end j farms
		
#----------------------------------------------------------------------------------------
# Solve all three solutions for each x
all.pred<-matrix();length(all.pred)<-3*n.sites*3*length(y)*2;dim(all.pred)<-c(3, n.sites, 3, length(y), 2)

for(j in 1:length(y)){# for each farm
	for(i in 1:n.sites){ #for each site
		for(s in 1:2){# for each species (Pink=1 and chum =2)
			
			#--------------------------------------
			#Copepodite
			#--------------------------------------
			# Left-hand solution
			all.pred[1,i,1,j,s]<-exp(b1*(x[i]-y[j]))*(1-exp(-b1*Lc[s]))*(a2+a1)/(b1*(a1-b1)*(a2+b1))-exp(a1*(x[i]-y[j]))*(1-exp(-a1*Lc[s]))*((b1+b2)/(a1*(a1+b2)*(a1-b1)))
	
			# Middle solution
			all.pred[1,i,2,j,s]<--(1-exp(a1*(x[i]-Lc[s]-y[j])))*(((b1+b2))/(a1*(a1+b2)*(a1-b1)))+(1-exp(b1*(x[i]-Lc[s]-y[j])))*((a2+a1)/(b1*(a1-b1)*(a2+b1)))+(1-exp(-b2*(x[i]-y[j])))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))-(1-exp(-a2*(x[i]-y[j])))*((b2+b1)/(a2*(a2+b1)*(a2-b2)))
			
			# Right hand solution
		all.pred[1,i,3,j,s]<-exp(-a2*(x[i]-y[j]))*(1-exp(a2*Lc[s]))*((b2+b1)/(a2*(a2-b2)*(a2+b1)))-exp(-b2*(x[i]-y[j]))*(1-exp(b2*Lc[s]))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))
			
			#--------------------------------------
			#Chalimus
			#--------------------------------------
			# Left-hand solution
			all.pred[2,i,1,j,s]<-exp(b1*(x[i]-y[j]))*(exp(-b1*Lc[s])-exp(-b1*Lc[s]*(1+Lh_R)))*(a2+a1)/(b1*(a1-b1)*(a2+b1)) -exp(a1*(x[i]-y[j]))*(exp(-a1*Lc[s])-exp(-a1*Lc[s]*(1+Lh_R)))*((b1+b2)/(a1*(a1+b2)*(a1-b1)))
	
			# Middle solution
			all.pred[2,i,2,j,s]<--(1-exp(a1*(x[i]-Lc[s]*(1+Lh_R)-y[j])))*(((b1+b2))/(a1*(a1+b2)*(a1-b1)))+(1-exp(b1*(x[i]-Lc[s]*(1+Lh_R)-y[j])))*((a2+a1)/(b1*(a1-b1)*(a2+b1)))+(1-exp(-b2*(x[i]-Lc[s]-y[j])))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))-(1-exp(-a2*(x[i]-Lc[s]-y[j])))* ((b2+b1)/(a2*(a2+b1)*(a2-b2)))
			
			# Right hand solution
		all.pred[2,i,3,j,s]<-exp(-a2*(x[i]-y[j]))*(exp(a2*Lc[s])-exp(a2*Lc[s]*(1+Lh_R)))*((b2+b1)/(a2*(a2-b2)*(a2+b1)))-exp(-b2*(x[i]-y[j]))*(exp(b2*Lc[s])-exp(b2*Lc[s]*(1+Lh_R)))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))
			
			#--------------------------------------
			#Motile
			#--------------------------------------
			# Left-hand solution
			all.pred[3,i,1,j,s]<-exp(b1*(x[i]-y[j]))*(exp(-b1*Lc[s]*(1+Lh_R))-exp(-b1*Lc[s]*(1+Lh_R+Lm_R)))*(a2+a1)/(b1*(a1-b1)*(a2+b1))-exp(a1*(x[i]-y[j]))*(exp(-a1*Lc[s]*(1+Lh_R))-exp(-a1*Lc[s]*(1+Lh_R+Lm_R)))*((b1+b2)/(a1*(a1+b2)*(a1-b1)))
	
			# Middle solution
			all.pred[3,i,2,j,s]<--(1-exp(a1*(x[i]-Lc[s]*(1+Lh_R+Lm_R)-y[j])))*(((b1+b2))/(a1*(a1+b2)*(a1-b1)))+(1-exp(b1*(x[i]-Lc[s]*(1+Lh_R+Lm_R)-y[j])))*((a2+a1)/(b1*(a1-b1)*(a2+b1)))+(1-exp(-b2*(x[i]-Lc[s]*(1+Lh_R)-y[j])))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))-(1-exp(-a2*(x[i]-Lc[s]*(1+Lh_R)-y[j])))*((b2+b1)/(a2*(a2+b1)*(a2-b2)))
			
			# Right hand solution
		all.pred[3,i,3,j,s]<-exp(-a2*(x[i]-y[j]))*(exp(a2*Lc[s]*(1+Lh_R))-exp(a2*Lc[s]*(1+Lh_R+Lm_R)))*((b2+b1)/(a2*(a2-b2)*(a2+b1)))-exp(-b2*(x[i]-y[j]))*(exp(b2*Lc[s]*(1+Lh_R))-exp(b2*Lc[s]*(1+Lh_R+Lm_R)))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))
	
		}# end species s
	} #end site i
} # end farm j 

#----------------------------------------------------------------------------------------
# Multiple by the index and sum all three solutions together
# should only be one solution that is non-zero for each location
int.solns<-matrix();length(int.solns)<-3*n.sites*length(y)*2; dim(int.solns)<-c(3,n.sites,length(y),2)
for(m in 1:3){ # for each stage
	for(i in 1:n.sites){ #for each site
		for(j in 1:length(y)){ #for each farm
			for(s in 1:2){ # for each species
					int.solns[m,i,j,s]<-(all.pred[m,i,1,j,s]*index[m,i,1,j,s] + all.pred[m,i,2,j,s]*index[m,i,2,j,s] + all.pred[m,i,3,j,s]*index[m,i,3,j,s])/c1
			}}}}
		
#----------------------------------------------------------------------------------------
# Add solutions from all three farms, scaled by their input
C_hat<-matrix(nrow=n.sites, ncol=2); H_hat<-C_hat; M_hat<-C_hat
for(i in 1:n.sites){ # for each site
	for(s in 1:2){ # for each species
		C_hat[i,s]<-k[s]*Lc[s] + phi[rep[i], 1]*int.solns[1,i,1,s] + phi[rep[i], 2]*int.solns[1,i,2,s] + phi[rep[i], 3]*int.solns[1,i,3,s]
		H_hat[i,s]<-sc[s]*(k[s]*Lc[s]*Lh_R + phi[rep[i], 1]*int.solns[1,i,1,s] + phi[rep[i], 2]*int.solns[1,i,2,s] + phi[rep[i], 3]*int.solns[1,i,3,s])
		M_hat[i,s]<-sh[s]*sc[s]*(k[s]*Lc[s]*Lm_R + phi[rep[i], 1]*int.solns[1,i,1,s] + phi[rep[i], 2]*int.solns[1,i,2,s] + phi[rep[i], 3]*int.solns[1,i,3,s])
	}}

# #-----------------------------------------------------------------------------------------------		
# #Confront with data
# C_prob<-r/(r+C_hat)
# H_prob<-r/(r+H_hat)
# M_prob<-r/(r+M_hat)
# for(i in 1:length(C)){
	# C[i] ~ dnegbin(C_prob[site[i],species[i]], r)
	# H[i] ~ dnegbin(H_prob[site[i],species[i]], r)
	# M[i] ~ dnegbin(M_prob[site[i],species[i]], r)
# } 

	return(list(C_hat, H_hat, M_hat))
	
	}