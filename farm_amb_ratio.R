library(here)
library(dclone)

load('~/Google Drive/Greens/FINALLY/GitHub/Spatiotemporal-infection-model/Workspaces/LiceFits2Data_20160215.RData')

S<-summary(mod)

kappa.log <- rnorm(1000, S[[1]]['k.log',1], S[[1]]['k.log',2])
par(mfrow=c(3,1))
#hist(kappa.log)
hist(exp(kappa.log))

phi.log <- rnorm(1000, S[[1]]['phi.log',1], S[[1]]['phi.log',2])
#hist(phi.log)
hist(exp(phi.log))

farmOverAmb <- exp(phi.log)/exp(kappa.log)
hist(farmOverAmb)
abline(v = exp(S[[1]]['phi.log',1])/exp(S[[1]]['k.log',1]), col=2)

quantile(farmOverAmb, c(0.025, 0.975))


############################

fPressure <- exp(S[[1]]['phi.log',1]) * farmL

range(fPressure / exp(S[[1]]['k.log',1]))

image(farmL)
