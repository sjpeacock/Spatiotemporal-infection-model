## This file is sourced by copes.R and simulates the larval distribution for one farm
## It uses multiple processors through the package snowfall

###############################################################################################
# 1) Define constants
###############################################################################################

#model parameters
D = 22.66667; 	# Diffusion
muN = 0.8;	 	# Naupliar removal (mortality or development)
muP = 0.2; 		# Copepodid removal
gamma = 1.56; 	# Advection

###############################################################################################
# 2) Nauplii distribution
###############################################################################################

#Simulate the convolution in time (space handled algebraically)
#for each timestep (starting at t=0 January 1, 2006), find the distribution of nauplii (area under curve = 1)

# Nauplii function for timestep j
Naup.fun<-function(j){
	Sum = 0;
	for (i in 0:(j-1)){
        Sum = Sum + Gn(x, (j-i)*dt, muN)*fset(Tsource,r1,r2,i*dt);
    	}#end i
	return(Sum*dt)}

#-----------------------------------------------------------------------------------------------
# Parallel computation
t.start<-proc.time()
sfInit(parallel=TRUE, cpus=cpus)
sfExport("muN", "dt", "x", "Tsource", "r1", "r2", "D", "gamma", "n", "f0")
sfSource("Larval-distribution/functions.R")
X<-sfLapply(1:n, Naup.fun)
sfStop()
cat("Nauplii process time (minutes) = ", (proc.time()-t.start)[3]/60)

Kn = matrix(0, nrow=length(x), ncol=n); #Kn is the array of naupliar density
for(i in 1:n){
	Kn[,i]<-X[[i]]}

#-----------------------------------------------------------------------------------------------
#scale distribution once decay starts to distribution at steady state
#Knpdf<-Kn/(sum(Kn[,zeroindex])*dx);


###############################################################################################
# 3) Copepidite distribution
###############################################################################################

#convolve the naupliar distribution with Green's function in space (using fourier transform) and time 

# Function to convolve naupliar distribution with Green's function
Cope.fun<-function(j){
	#integrate the convolution from the startpoint to iteration*dt in t to get the copepodid density
	Sum = 0;
	for (i in 1:j){
		Kx = Re(fft(fft(t(Gn(x, (j-(i-1))*dt, muP)))*fft(Kn[,i]), inverse=TRUE)); #the convolution
		Kx = c(Kx[(which(x==0)+1):(length(x))], Kx[1:which(x==0)])*dx; #account for frustrating 
	    Sum = Sum + Kx;
	   }#end i
	return(Sum*dt);
	}
	
#-----------------------------------------------------------------------------------------------

# Parallel computation for copepodid distribution (takes 30 minutes)
t.start<-proc.time()
sfInit(parallel=TRUE, cpus=cpus)
sfExport("D","gamma","muN", "muP", "x", "Kn", "dx", "dt", "Tsource", "r1", "r2", "n", "f0")
sfSource("Larval-distribution/functions.R")
X<-sfLapply(x=c(1:n), fun=Cope.fun)
sfStop()
cat("Process time (minutes) = ", (proc.time()-t.start)[3]/60)

K = matrix(0, nrow=length(x), ncol=n); #Kn is the array of naupliar density
for(i in 1:n) K[,i]<-X[[i]]

#-----------------------------------------------------------------------------------------------

#Kpdf = K/(sum(K[,zeroindex])*dx);
    

#------------------------------------------------------------------------------------------------------------

est = matrix(0, nrow=length(x), ncol=n);

if (Xsource>0){
    est[(round(Xsource/dx)+1):length(x),] = K[1:(length(x)-round(Xsource/dx)),]}else if (Xsource<0){
    est[1:(length(x)-round(-Xsource/dx)),] = K[(round(-Xsource/dx)+1):length(x),]} else {
    est = K}

