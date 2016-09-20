simulate.lice<-function(param, dist, day, log=TRUE){
	
	#Fixed parameters
	tauC<-0.525^-2
	tauH<-0.250^-2
	tauM<-0.187^-2
	
	#Indices for summing over time (do not depend on free parameters)
	C_index<-c(0:(round(tauC/dt)-1))
	H_index<-c(round(tauC/dt):(round(tauH/dt)-1))
	M_index<-c(round(tauH/dt):(round(tauM/dt)-1))
	

n.sites<-length(dist)
	
	if(log==TRUE){
		k.log<-param['k.log']; k <- exp(k.log)
		phi.log<-param['phi.log']; phi<-exp(phi.log)
		V.log<-param['V.log'];V <- exp(V.log)
		Psc.log<-param['Psc.log']; Psc<- exp(Psc.log)/(1+exp(Psc.log))
		Psh.log<-param['Psh.log']; Psh<- exp(Psh.log)/(1+exp(Psh.log))
		Csc.log<-param['Csc.log']; Csc<- exp(Csc.log)/(1+exp(Csc.log))
		Csh.log<-param['Csh.log'];Csh<- exp(Csh.log)/(1+exp(Csh.log))
		r.log<-param['r.log'];r<-exp(r.log)
	}else{
		k<-param['k']
		phi<-param['phi']
		V<-param['V']
		Psc<-param['Psc']
		Psh<-param['Psh']
		Csc<-param['Csc']
		Csh<-param['Csh']
		r<-param['r']
	}	
	
		
	Vind <- V*dt/dx;
    t.C<-matrix(nrow=n.sites, ncol=length(C_index)); x.C<-t.C
    t.H<-matrix(nrow=n.sites, ncol=length(H_index)); x.H<-t.H
    t.M<-matrix(nrow=n.sites, ncol=length(M_index)); x.M<-t.M
    for(i in 1:n.sites){
    	for(j in 1:length(C_index)){ 
    		t.C[i,j]<-round(day[i]/dt+t0-C_index[j])
    		x.C[i,j]<-max(0, round(dist[i]/dx+x0-Vind*C_index[j]))+1} #!! Add a row of zeros to the x side of L
    	for(j in 1:length(H_index)){
    		t.H[i,j]<-round(day[i]/dt+t0-H_index[j])
    		x.H[i,j]<-max(0, round(dist[i]/dx+x0-Vind*H_index[j]))+1}
    	for(j in 1:length(M_index)){
    		t.M[i,j]<-round(day[i]/dt+t0-M_index[j])
    		x.M[i,j]<-max(0, round(dist[i]/dx+x0-Vind*M_index[j]))+1}
    	}
    C_int<-matrix(nrow=n.sites, ncol=length(C_index)); C_hat<-matrix(nrow=n.sites, ncol=2)
    H_int<-matrix(nrow=n.sites, ncol=length(H_index)); H_hat<-matrix(nrow=n.sites, ncol=2)
    M_int<-matrix(nrow=n.sites, ncol=length(M_index)); M_hat<-matrix(nrow=n.sites, ncol=2)
    for(i in 1:n.sites){ 
    	for(j in 1:length(C_index)){
        		C_int[i,j] <- L[x.C[i,j],t.C[i,j]] }
        	C_hat[i,1] <- k*tauC+phi*sum(C_int[i,])*dt;
			C_hat[i,2] <- k*tauC+phi*sum(C_int[i,])*dt;
			for(j in 1:length(H_index)){
				H_int[i,j] <- L[x.H[i,j],t.H[i,j]] }
			H_hat[i,1]<- Psc*(k*tauH+phi*sum(H_int[i,])*dt);
			H_hat[i,2]<- Csc*(k*tauH+phi*sum(H_int[i,])*dt);
			for(j in 1:length(M_index)){
				M_int[i,j] <- L[x.M[i,j],t.M[i,j]] }
			M_hat[i,1] <- Psc*Psh*(k*tauM+phi*sum(M_int[i,])*dt);
			M_hat[i,2] <- Csc*Csh*(k*tauM+phi*sum(M_int[i,])*dt);		
		} 
	return(list(C_hat, H_hat, M_hat))
}

#########################################################
logLike.nb<-function(param){
	
	pred<-simulate.lice(param, dat$dist, dat$day)
	 
	r<-exp(param['r.log'])
	 	
  	LL<-sum(dnbinom(dat$C, mu=pred[[1]][dat$site.id,dat$spp[i]], size=r, log=TRUE),dnbinom(dat$H, mu=pred[[2]][dat$site.id,dat$spp[i]], size=r, log=TRUE),dnbinom(dat$M, mu=pred[[3]][dat$site.id,dat$spp[i]], size=r, log=TRUE))
  		
  	
  	return(-LL)	
  	}
  	 

