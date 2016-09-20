###############
# Copedist Domain
################
#-----------------------------------------------------------------------------------------------
#space and time domains
dx = 0.05							# resolution of spatial grid
Xmax = 80							# maximum X to the left of Shewell Island
x = c((-60/dx):(80/dx-dx))*dx	# spatial grid

#-----------------------------------------------------------------------------------------------
#time positions 
dt = 0.05							# resolution of temporal grid
Tmax = 150
n = round(Tmax/dt)					#number of grid spaces since SLICE treatment
T = c((10/dt):n)*dt

x0=which(x==0);
t0=(-10/dt)
