####Script to Create SIBER Ellipses####

##Packages
if (!require("SIBER")) install.packages("SIBER")
library(SIBER) #Required for generating physiological ellipses and percent overlap calculations

#Note: Also requires JAGS downloaded on local computer

##Set working directory
getwd()
setwd("C:/Users/snhac/OneDrive/Desktop/Bonaire_memory/Outputs/SIBER")

##Read in SIBER Object 
#This was created in 05_PhysiologySIBER.Rmd
load("Phys.SIBER.Object.RData")

##set up the model parameters
#options for running jags
parms <- list()
parms$n.iter <- 2 * 10^6   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^4 # discard the first set of values
parms$n.thin <- 100     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
parms$save.output = T
parms$save.dir = getwd()

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3


# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the
# means. Fitting is via the JAGS method.
phys.ellipses.posterior <- siberMVN(phys.siber.object, parms, priors)

#calculate posterior estimates of ellipses for all groups
SEA.B<-siberEllipses(phys.ellipses.posterior)
#set column names for SEA.B
colnames(SEA.B)<-names(phys.ellipses.posterior)

##Save Ellipse Data in Outputs
save(phys.ellipses.posterior, SEA.B, file="Phys.SIBER.Ellipses.RData")
