###############
# Copedist Domain
################
#-----------------------------------------------------------------------------------------------
#space and time domains
dx = 0.5							# resolution of spatial grid
Xmax = 80							# maximum X to the left of Shewell Island
x = c((-60/dx):(80/dx-dx))*dx	# spatial grid

#-----------------------------------------------------------------------------------------------
#time positions 
dt = 0.5							# resolution of temporal grid
Tmax = 303							# 303 is July 1, 2006, assuming t=0 is Sept 1, 2005
n = round(Tmax/dt)					#number of grid spaces since SLICE treatment
T = c(0:(n-1))*dt

x0=which(x==0);
t0=which(T==0)
