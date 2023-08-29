#Title: "Calculation of Physiological and Performance Metrics"
#Author: "Serena Hackerott"
#Date: "08/29/2023"

#-------Set Up------------

####Load Packages####

#Note: Run "Graphing Parameters" section from 01_ExperimentalSetup.R file


####Load Data####

Prot<-read.csv("Data/Protein.csv", header=TRUE)
Bio<-read.csv("Data/Biomass.csv", header=TRUE)

Sym<-read.csv("Data/Symbionts.csv", header=TRUE)
Chl<-read.csv("Data/Chlorophyll.csv", header=TRUE)
PAM<-read.csv("Data/PAMData.csv", header=TRUE)
PAM_Meta<-read.csv("Data/PAMMeta.csv", header=TRUE)

BCA.Stand<-read.csv("Data/BCAStandards.csv", header=TRUE)
Wax.Stand<-read.csv("Data/WaxStandards.csv", header=TRUE)

SampData<-read.csv("Data/Samples.csv", header=TRUE)


####Sample Meta Data####
str(SampData)



##Set factor variables

##Add a Sample Set Variable

####Sample Surface Area####

##Calculate Surface Area of Cylinders from Standards
Wax.Stand$SA_cm2<-2*pi*(Wax.Stand$Dm_cm/2)*Wax.Stand$Ht_cm+2*pi*(Wax.Stand$Dm_cm/2)^2

##Calculate Difference in Wax Weight
Wax.Stand$Wax.D_g<- Wax.Stand$Wax.F_g-Wax.Stand$Wax.I_g

##Linear Model of SA as a Function of Wax Weight
wax.SA.lm  <- lm(SA_cm2~Wax.D_g, data=Wax.Stand)
summary(wax.SA.lm)
coef(wax.SA.lm)

SA.mod <- function(wax.weight) {
  coefs <- coef(wax.SA.lm)
  #y = mx + b
  SA <- ((coefs[2] * wax.weight) + coefs[1])
  return(SA)}

##Calculate Difference in Wax Weight (g) in Sample Data
Samp$Wax.D_g<- Samp$Wax.F_g-Samp$Wax.I_g

##Calculate Surface Area (cm^2) by Applying Linear Model Function
Samp$SA_cm2<-SA.mod(Samp$Wax.D_g)


#-----Physiology Metrics------------

