#########################################################
model<-function(){
	Vind <- V*dt/dx;
    #t.C<-matrix(nrow=n.sites, ncol=length(C_index)); x.C<-t.C
    #t.H<-matrix(nrow=n.sites, ncol=length(H_index)); x.H<-t.H
    #t.M<-matrix(nrow=n.sites, ncol=length(M_index)); x.M<-t.M
    
    
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
    
    # C_int<-matrix(nrow=n.sites, ncol=length(C_index)); C_hat<-matrix(nrow=n.sites, ncol=2)
    # H_int<-matrix(nrow=n.sites, ncol=length(H_index)); H_hat<-matrix(nrow=n.sites, ncol=2)
    # M_int<-matrix(nrow=n.sites, ncol=length(M_index)); M_hat<-matrix(nrow=n.sites, ncol=2)
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
  	
  	#par(mfrow=c(3,1)); plot(dist, C_hat[,1]); plot(dist, H_hat[,1]); plot(dist, M_hat[,1])
  	C_prob<-r/(r+C_hat)
  	H_prob<-r/(r+H_hat)
  	M_prob<-r/(r+M_hat)
  	for(i in 1:length(C)){
  		C[i] ~ dnegbin(C_prob[site.id[i],spp[i]], r)
  		H[i] ~ dnegbin(H_prob[site.id[i],spp[i]], r)
  		M[i] ~ dnegbin(M_prob[site.id[i],spp[i]], r)
  		} 
 	k.log ~ dnorm(mu[1], tau[1]);k <- exp(k.log)
	phi.log ~ dnorm(mu[2], tau[2]);phi<-exp(phi.log)
	V.log ~ dnorm(mu[3], tau[3]); V <- exp(V.log)	
	Psc.log ~ dnorm(mu[4], tau[4]); Psc<- exp(Psc.log)/(1+exp(Psc.log))
	Psh.log ~ dnorm(mu[5], tau[5]);	Psh<- exp(Psh.log)/(1+exp(Psh.log))
	Csc.log ~ dnorm(mu[6], tau[6]);Csc<- exp(Csc.log)/(1+exp(Csc.log))
	Csh.log ~ dnorm(mu[7], tau[7]);Csh<- exp(Csh.log)/(1+exp(Csh.log))
	r.log ~ dnorm(mu[8], tau[8]); r<-exp(r.log)
}

