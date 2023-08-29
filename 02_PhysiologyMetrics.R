#Title: "Calculation of Physiological and Performance Metrics"
#Author: "Serena Hackerott"
#Date: "08/29/2023"

#-------Set Up------------

####Load Packages####

library(plyr)
library(ggplot2)

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
SampData$TimeP<-factor(SampData$TimeP, levels=c("W2", "M1", "M4", "M8", "M12"), ordered=TRUE)
SampData$Site<-factor(SampData$Site, levels=c("SS", "KL"), ordered=TRUE)
SampData$Genotype<-factor(SampData$Genotype, levels=c("AC8", "AC10", "AC12"), ordered=TRUE)
SampData$Orig<-factor(SampData$Orig, levels=c("N", "T"), ordered=TRUE)
SampData$Origin<-factor(SampData$Origin, levels=c("Native", "Transplant"), ordered=TRUE)

##Add a Sample Set Variable
SampData$Set<-paste(SampData$TimeP, SampData$Site, SampData$Genotype, SampData$Treat, SampData$Orig, sep=".")



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
SampData$Wax.D_g<- SampData$Wax.F_g-SampData$Wax.I_g

##Calculate Surface Area (cm^2) by Applying Linear Model Function
SampData$SA_cm2<-SA.mod(SampData$Wax.D_g)


#-----Physiology Metrics------------

####Subset Sample MetaData by Experiment
Corals_Main<-subset(SampData, Treat=="M")
Corals_Therm<-subset(SampData, Treat=="C" | Treat=="H")


####Protein####
str(Prot)

####Standard Curve

##Subset Standards Data
Prot.St<- subset(Prot, Input=="Standard")

##Merge Standard Meta Data with Absorbance data
#Merges by Random Number (RandN) column
BCAStands<-merge(BCA.Stand, Prot.St, all.x=TRUE)

##Plot Protein Concentration as a function of Absorbance 
plot(BCAStands$A562, BCAStands$Protein_ug.ml, pch=19)

##Fit Models of Protein Concentration as a function of Absorbance 

#First degree polynomial equation:
BCA.lm  <- lm(Protein_ug.ml~A562, data=BCAStands)
#Second degree polynomial equation:
BCA.ploy2 <- lm(Protein_ug.ml~poly(A562,2,raw=TRUE), data=BCAStands)
#Third degree polynomial equation:
BCA.ploy3 <- lm(Protein_ug.ml~poly(A562,3,raw=TRUE), data=BCAStands)

##Add Models to Plot
xx <- seq(0,3, length=50)

lines(xx, predict(BCA.lm, newdata=data.frame(A562=xx)), col="red")
lines(xx, predict(BCA.ploy2, newdata=data.frame(A562=xx)), col="green")
lines(xx, predict(BCA.ploy3, newdata=data.frame(A562=xx)), col="blue")

##Compare Models for Best Fit
anova(BCA.lm, BCA.ploy2)
#Poly2 is a significantly better fit than linear

anova(BCA.ploy2, BCA.ploy3)
#Poly3 is significantly better fit than Ploy2

##Create a function with the Poly3 model Equation
summary(BCA.ploy3)
coef(BCA.ploy3)

TP.mod <- function(absorbance) {
  coefs <- coef(BCA.ploy3)
  #y = d + cx + bx^2 + ax^3
  protein <- coefs[1] + (coefs[2] * absorbance) + (coefs[3] * absorbance^2) + (coefs[4] * absorbance^3)
  return(protein)}


####Calculate Protein Concentration

##Subset Sample Data (Unknown Sample)
Prot.Un<- subset(Prot, Input=="Un")

##Calculate Protein Concentration (ug/ml)
Prot.Un$TP_ug.ml<-TP.mod(Prot.Un$A562)

##Merge with Sample Meta Data to Calculate Protein per Surface Area
#Merges by Random Number (RandN) column
#Retains Sample Meta Data only for samples with Protein data (Main samples only)
Prot.Un<-merge(Prot.Un, SampData, all.x=TRUE, all.y=FALSE)

##Calculate Total Protein (ug) 
Prot.Un$TP_ug<-Prot.Un$TP_ug.ml*Prot.Un$Vol_ml

##Calculate Protein per Surface Area (ug/cm^2)
Prot.Un$TP_ug.cm2<-Prot.Un$TP_ug/Prot.Un$SA_cm2

##Separate Coral and Symbiont Fractions
Prot.C<- subset(Prot.Un, Fraction=="C")
Prot.C<-rename(Prot.C, c("TP_ug.cm2" = "TP_ug.cm2_C"))

Prot.S<- subset(Prot.Un, Fraction=="S")
Prot.S<-rename(Prot.S, c("TP_ug.cm2" = "TP_ug.cm2_S"))

####Check for Outliers

##Check for Outliers in Protein of Coral Host by Timepoint

#Week 2
ggplot(subset(Prot.C, TimeP=="W2"), aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1000)+
  theme(axis.text.x = element_text(angle = 90))

#Month 1
ggplot(subset(Prot.C, TimeP=="M1"), aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1000)+
  theme(axis.text.x = element_text(angle = 90))

Prot.C$RandN[which(Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C>600)] #"15" "15" "15" "29" "29" "29"
Prot.C$RandN[which(Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C>500 & Prot.C$Site=="KL")] #"4"
Prot.C$RandN[which(Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C<200)] #"20" "20" "20"

#Month 4
ggplot(subset(Prot.C, TimeP=="M4"), aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1000)+
  theme(axis.text.x = element_text(angle = 90))

Prot.C$RandN[which(Prot.C$TimeP=="M4" & Prot.C$TP_ug.cm2_C<170)] #"121" "125" "126" "127"
Prot.C$RandN[which(Prot.C$TimeP=="M4" & Prot.C$TP_ug.cm2_C<250 & Prot.C$Site=="SS")] #"125" "126" "140" "140"

#Month 8
ggplot(subset(Prot.C, TimeP=="M8"), aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(100,1000)+
  theme(axis.text.x = element_text(angle = 90))

Prot.C$RandN[which(Prot.C$TimeP=="M8" & Prot.C$TP_ug.cm2_C>550)] #"184" "184" "186" "187" "188"

#Month 12
ggplot(subset(Prot.C, TimeP=="M12"), aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1100)+
  theme(axis.text.x = element_text(angle = 90))

Prot.C$RandN[which(Prot.C$TimeP=="M12" & Prot.C$TP_ug.cm2_C>600)] #"211" "222" "222" "222"

##Remove Outlier Readings
Prot.C.o<-Prot.C[-c(which((Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C>600) |
(Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C>500 & Prot.C$Site=="KL") |
(Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C<200) |
(Prot.C$TimeP=="M4" & Prot.C$TP_ug.cm2_C<170) |
(Prot.C$TimeP=="M4" & Prot.C$TP_ug.cm2_C<250 & Prot.C$Site=="SS") |
(Prot.C$TimeP=="M8" & Prot.C$TP_ug.cm2_C>550) |
(Prot.C$TimeP=="M12" & Prot.C$TP_ug.cm2_C>600))),]

ggplot(Prot.C.o, aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(100,600)+
  theme(axis.text.x = element_text(angle = 90))


##Check for Outliers in Protein of Symbiont by Timepoint

#Week 2
ggplot(subset(Prot.S, TimeP=="W2"), aes(x=Set, y=TP_ug.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50, 1200)+
  theme(axis.text.x = element_text(angle = 90))

Prot.S$RandN[which(Prot.S$TimeP=="W2" & Prot.S$TP_ug.cm2_S>750)] #"49" "52" "84"

#Month 1
ggplot(subset(Prot.S, TimeP=="M1"), aes(x=Set, y=TP_ug.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1400)+
  theme(axis.text.x = element_text(angle = 90))

Prot.S$RandN[which(Prot.S$TimeP=="M1" & Prot.S$TP_ug.cm2_S>750)] #"15" "4"  "41" "47" "5"
Prot.S$RandN[which(Prot.S$TimeP=="M1" & Prot.S$TP_ug.cm2_S<250)] #"20" "20" "20"

#Month 4
ggplot(subset(Prot.S, TimeP=="M4"), aes(x=Set, y=TP_ug.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1000)+
  theme(axis.text.x = element_text(angle = 90))

Prot.S$RandN[which(Prot.S$TimeP=="M4" & Prot.S$TP_ug.cm2_S>500 & Prot.S$Site=="SS" & Prot.S$Genotype=="AC8")] #"129"

#Month 8
ggplot(subset(Prot.S, TimeP=="M8"), aes(x=Set, y=TP_ug.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1000)+
  theme(axis.text.x = element_text(angle = 90))

#Month 12
ggplot(subset(Prot.S, TimeP=="M12"), aes(x=Set, y=TP_ug.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1000)+
  theme(axis.text.x = element_text(angle = 90))

Prot.S$RandN[which(Prot.S$TimeP=="M12" & Prot.S$TP_ug.cm2_S>750)] #"222" "222" "222"

##Remove Outlier Readings
Prot.S.o<-Prot.S[-c(which((Prot.S$TimeP=="W2" & Prot.S$TP_ug.cm2_S>750) |
(Prot.S$TimeP=="M1" & Prot.S$TP_ug.cm2_S>750) |
(Prot.S$TimeP=="M1" & Prot.S$TP_ug.cm2_S<250) |
(Prot.S$TimeP=="M4" & Prot.S$TP_ug.cm2_S>500 & Prot.S$Site=="SS" & Prot.S$Genotype=="AC8") |
(Prot.S$TimeP=="M12" & Prot.S$TP_ug.cm2_S>750))),]

ggplot(Prot.S.o, aes(x=Set, y=TP_ug.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(100,750)+
  theme(axis.text.x = element_text(angle = 90))


####Average Across Replicate Readings ("Rep" column: A, B, C)
names(Prot.C.o)
Prot.C.a<-aggregate(Prot.C.o$TP_ug.cm2_C, list(Prot.C.o$RandN, Prot.C.o$ID), mean)
names(Prot.C.a)<-c("RandN", "ID", "TP_ug.cm2_C")

names(Prot.S.o)
Prot.S.a<-aggregate(Prot.S.o$TP_ug.cm2_S, list(Prot.S.o$RandN, Prot.S.o$ID), mean)
names(Prot.S.a)<-c("RandN", "ID", "TP_ug.cm2_S")


####Add Total Protein to Main Coral Data

##Merge Averaged Protein Data of Coral Host (C) and Symbiont (S) with Coral Main Data
#Merges by Random Number (RandN) column
#Retains all Main samples
#Retains RandN, ID, TimeP, Site, Genotype, Orig, Origin, Set, and SA_cm2 from Corals_Main
names(Corals_Main)

CoralData<-merge(Corals_Main[,c(1:5, 8:9, 13, 15)], Prot.C.a, all.x=TRUE)
CoralData<-merge(CoralData, Prot.S.a, all.x=TRUE)


####Biomass####


####Symbionts####


####Chlorophyll####


####Fv/Fm####

