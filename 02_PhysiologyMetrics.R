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
#Adds necessary Slurry Volume (Vol_ml) and Surface Area (SA_cm2) columns
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

##Initial Visual Check
ggplot(Prot.C, aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))

ggplot(Prot.S, aes(x=Set, y=TP_ug.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))


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

Prot.C$RandN[which(Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C<200)] #"20" "20" "20"

#Month 4
ggplot(subset(Prot.C, TimeP=="M4"), aes(x=Set, y=TP_ug.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(50,1000)+
  theme(axis.text.x = element_text(angle = 90))

Prot.C$RandN[which(Prot.C$TimeP=="M4" & Prot.C$TP_ug.cm2_C<170)] #"121" "125" "126" "127"

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
Prot.C.o<-Prot.C[-c(which((Prot.C$TimeP=="M1" & Prot.C$TP_ug.cm2_C<200) |
(Prot.C$TimeP=="M4" & Prot.C$TP_ug.cm2_C<170) |
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

Prot.S$RandN[which(Prot.S$TimeP=="M1" & Prot.S$TP_ug.cm2_S>1000)] #"41"
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
(Prot.S$TimeP=="M1" & Prot.S$TP_ug.cm2_S>1000) |
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
#Merges by Random Number (RandN) and ID columns
#Retains all Main samples
#Retains RandN, ID, TimeP, Site, Genotype, Orig, Origin, Set, and SA_cm2 from Corals_Main
names(Corals_Main)

CoralData<-merge(Corals_Main[,c(1:5, 8:9, 13, 15)], Prot.C.a, all.x=TRUE)
CoralData<-merge(CoralData, Prot.S.a, all.x=TRUE)


####Biomass####
str(Bio)

####Calculate Ash Free Dry Weight (AFDW)

##Calculate Change in Weight after Burning (g)
Bio$DeltaBurn_g<-Bio$Dried_g-Bio$Burned_g

##Calculate Ash Free Dry Weight (AFDW) (g) per Volume Input of Tissue (ml) 
Bio$AFDW_g.ml<-Bio$DeltaBurn_g/Bio$InVol_ml

##Merge with Sample Meta Data to Calculate Biomass per Surface Area
#Merges by Random Number (RandN) column
#Adds necessary Slurry Volume (Vol_ml) and Surface Area (SA_cm2) columns
#Retains Sample Meta Data only for samples with Biomass data (Main samples only)
Bio<-merge(Bio, SampData, all.x=TRUE, all.y=FALSE)

##Calculate Total AFDW (g) 
Bio$AFDW_g<-Bio$AFDW_g.ml*Bio$Vol_ml

##Calculate AFDW per Surface Area (mg/cm^2)
Bio$AFDW_g.cm2<-Bio$AFDW_g/Bio$SA_cm2
Bio$AFDW_mg.cm2 <- Bio$AFDW_g.cm2 * 1000

##Separate Coral and Symbiont Fractions
Bio.C<-subset(Bio, Fraction=="C")
Bio.C<-rename(Bio.C, c("AFDW_mg.cm2" = "AFDW_mg.cm2_C"))

Bio.S<-subset(Bio, Fraction=="Z")
Bio.S<-rename(Bio.S, c("AFDW_mg.cm2" = "AFDW_mg.cm2_S"))


##Initial Visual Check
ggplot(Bio.C, aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))

ggplot(Bio.S, aes(x=Set, y=AFDW_mg.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))


####Check for Outliers

##Check for Outliers in Biomass of Coral Host by Timepoint

#Week 2
ggplot(subset(Bio.C, TimeP=="W2"), aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 1
ggplot(subset(Bio.C, TimeP=="M1"), aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))

Bio.C$RandN[which(Bio.C$TimeP=="M1" & Bio.C$AFDW_mg.cm2_C<0.55)] #20 #Also outlier for Protein

#Month 4
ggplot(subset(Bio.C, TimeP=="M4"), aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))

Bio.C$RandN[which(Bio.C$TimeP=="M4" & Bio.C$AFDW_mg.cm2_C>20)] #123 #Need to check this

ggplot(subset(Bio.C, TimeP=="M4"), aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 8
ggplot(subset(Bio.C, TimeP=="M8"), aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 12
ggplot(subset(Bio.C, TimeP=="M12"), aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,3.5)+
  theme(axis.text.x = element_text(angle = 90))

Bio.C$RandN[which(Bio.C$TimeP=="M12" & Bio.C$AFDW_mg.cm2_C>3)] #222 #Also outlier for Protein

##Remove Outliers 
Bio.C.o<-Bio.C[-c(which((Bio.C$TimeP=="M1" & Bio.C$AFDW_mg.cm2_C<0.55) |
(Bio.C$TimeP=="M4" & Bio.C$AFDW_mg.cm2_C>20) |
(Bio.C$TimeP=="M12" & Bio.C$AFDW_mg.cm2_C>3))),]

ggplot(Bio.C.o, aes(x=Set, y=AFDW_mg.cm2_C)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))


##Check for Outliers in Biomass of Symbiont by Timepoint

#Week 2
ggplot(subset(Bio.S, TimeP=="W2"), aes(x=Set, y=AFDW_mg.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 1
ggplot(subset(Bio.S, TimeP=="M1"), aes(x=Set, y=AFDW_mg.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(-0.5,1.5)+
  theme(axis.text.x = element_text(angle = 90))

Bio.S$RandN[which(Bio.S$TimeP=="M1" & Bio.S$AFDW_mg.cm2_S<0.05)] #25

#Month 4
ggplot(subset(Bio.S, TimeP=="M4"), aes(x=Set, y=AFDW_mg.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 8
ggplot(subset(Bio.S, TimeP=="M8"), aes(x=Set, y=AFDW_mg.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 12
ggplot(subset(Bio.S, TimeP=="M12"), aes(x=Set, y=AFDW_mg.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 90))

Bio.S$RandN[which(Bio.S$TimeP=="M12" & Bio.S$AFDW_mg.cm2_S>1.5)] #222 #Also outlier for Protein and Biomass of Host


##Remove Outliers 
Bio.S.o<-Bio.S[-c(which((Bio.S$TimeP=="M1" & Bio.S$AFDW_mg.cm2_S<0.05) |
  (Bio.S$TimeP=="M12" & Bio.S$AFDW_mg.cm2_S>1.5))),]

ggplot(Bio.S.o, aes(x=Set, y=AFDW_mg.cm2_S)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1.5)+
  theme(axis.text.x = element_text(angle = 90))


####Add AFDW of Coral and Symbiont to Main Coral Data

##Merge Cleaned Biomass Data of Coral Host (C) and Symbiont (S) with Coral Main Data
#Merges by ID column and adds AFDW (mg/cm^2) (AFDW_mg.cm2_C or _S) columns from Biomass dataframes
#Retains all Main samples
names(Bio.C.o)
CoralData<-merge(CoralData, Bio.C.o[,c(1, 25)], all.x=TRUE)
names(Bio.S.o)
CoralData<-merge(CoralData, Bio.S.o[,c(1, 25)], all.x=TRUE)

##Calculate Symbiont:Host Ratio of Biomass
CoralData$AFDW_mg.cm2_S.C<-CoralData$AFDW_mg.cm2_S/CoralData$AFDW_mg.cm2_C


####Chlorophyll####
str(Chl)

####Calculate Chlorophyll Concentration

##Equations for Dinos from Jeffrey and Humphrey 1975 in 100% acetone
#chla = 11.43*A663 - 0.64*A630
#chlc2 = 27.09*A630 - 3.63*A663

##Subtract Background A750 from A630 and A663
Chl$A630.c<-Chl$A630-Chl$A750
Chl$A663.c<-Chl$A663-Chl$A750

##Divide by Pathlength (0.5cm pathlength for 175ul sample in UVStar Plate) 
Chl$A630.c<-c(Chl$A630.c/0.5)
Chl$A663.c<-c(Chl$A663.c/0.5)

##Calculate Chl-a and Chl-c2 in µg/ml
Chl$Chl.a_ug.ml<-11.43*Chl$A663.c- 0.64*Chl$A630.c
Chl$Chl.c2_ug.ml <- 27.09*Chl$A630.c - 3.63*Chl$A663.c

##Merge with Sample Data to Calculate Chlorophyll per Surface Area
#Merges by Random Number (RandN) column
#Adds necessary Slurry Volume (Vol_ml) and Surface Area (SA_cm2) columns
Chl<-merge(Chl, SampData,  all=TRUE)

##Calculate Total Chlorophyll-a and c2 (ug) 
Chl$Chl.a_ug<-Chl$Chl.a_ug.ml*Chl$Vol_ml
Chl$Chl.c2_ug<-Chl$Chl.c2_ug.ml*Chl$Vol_ml

##Calculate Chlorophyll-a and c2 per Surface Area (ug/cm^2)
Chl$Chl.a_ug.cm2<-Chl$Chl.a_ug/Chl$SA_cm2
Chl$Chl.c2_ug.cm2<-Chl$Chl.c2_ug/Chl$SA_cm2

##Initial Visual Check
ggplot(Chl, aes(x=Set, y=Chl.a_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))

ggplot(Chl, aes(x=Set, y=Chl.c2_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))

##Set Negative Values to Zero
Chl$Chl.a_ug.cm2[(which(Chl$Chl.a_ug.cm2<0))]<-0
Chl$Chl.c2_ug.cm2[(which(Chl$Chl.c2_ug.cm2<0))]<-0

##Add Chlorophyll-a and c2 for Total Chlorophyll per Surface Area (ug/cm^2)
Chl$Chl_ug.cm2<-Chl$Chl.a_ug.cm2+Chl$Chl.c2_ug.cm2
  
ggplot(Chl, aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))


####Check for Outliers

##Check for Outliers in Chlorphyll by Treatment and Timepoint

##Main
Chl.M<-subset(Chl, Treat=="M")

#Week 2
ggplot(subset(Chl.M, TimeP=="W2"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 1
ggplot(subset(Chl.M, TimeP=="M1"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,3.5)+
  theme(axis.text.x = element_text(angle = 90))

Chl.M$RandN[which(Chl.M$TimeP=="M1" & Chl.M$Chl_ug.cm2<0.5)] #"20" "20" "20" #Also outlier for Protein and Biomass

#Month 4
ggplot(subset(Chl.M, TimeP=="M4"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,4.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 8
ggplot(subset(Chl.M, TimeP=="M8"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,4.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 12
ggplot(subset(Chl.M, TimeP=="M12"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,7.5)+
  theme(axis.text.x = element_text(angle = 90))

Chl.M$RandN[which(Chl.M$TimeP=="M12" & Chl.M$Chl_ug.cm2>6)] #"222" "222" "222" #Also outlier for Protein and Biomass


##Control
Chl.C<-subset(Chl, Treat=="C")

#Week 2
ggplot(subset(Chl.C, TimeP=="W2"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(-0.5,10)+
  theme(axis.text.x = element_text(angle = 90))

Chl.C$RandN[which(Chl.C$TimeP=="W2" & Chl.C$Chl_ug.cm2>2.5)] #"W2_24" "W2_8" 
##Need to Check Week 2

#Month 1
ggplot(subset(Chl.C, TimeP=="M1"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 4
ggplot(subset(Chl.C, TimeP=="M4"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.5)+
  theme(axis.text.x = element_text(angle = 90))

#Month 8
ggplot(subset(Chl.C, TimeP=="M8"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,3.5)+
  theme(axis.text.x = element_text(angle = 90))

Chl.C$RandN[which(Chl.C$TimeP=="M8" & Chl.C$Chl_ug.cm2>2.5)] #"M8_33" "M8_33" "M8_33" 
Chl.C$RandN[which(Chl.C$TimeP=="M8" & Chl.C$Chl_ug.cm2<0.5)] #"M8_12" "M8_12" "M8_12"

#Month 12
ggplot(subset(Chl.C, TimeP=="M12"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,4)+
  theme(axis.text.x = element_text(angle = 90))


##Heated
Chl.H<-subset(Chl, Treat=="H")
Chl.H<-Chl.H[-c(which(is.na(Chl.H$Chl_ug.cm2))),]

#Week 2
ggplot(subset(Chl.H, TimeP=="W2"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90))

Chl.H$RandN[which(Chl.H$TimeP=="W2" & Chl.H$Chl_ug.cm2>0.75)] #"W2_88" "W2_88" "W2_88"

#Month 1
ggplot(subset(Chl.H, TimeP=="M1"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90))

#Month 4
ggplot(subset(Chl.H, TimeP=="M4"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90))

#Month 8
ggplot(subset(Chl.H, TimeP=="M8"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 90))

Chl.H$RandN[which(Chl.H$TimeP=="M8" & Chl.H$Chl_ug.cm2>1.5)] #"M8_11" "M8_11" "M8_42" "M8_42"

#Month 12
ggplot(subset(Chl.H, TimeP=="M12"), aes(x=Set, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90))

Chl.H$RandN[which(Chl.H$TimeP=="M12" & Chl.H$Chl_ug.cm2>0.51)] #M12_63" "M12_63"


##Remove Outliers 
Chl.M.o<-Chl.M[-c(which((Chl.M$TimeP=="M1" & Chl.M$Chl_ug.cm2<0.5) | 
                          (Chl.M$TimeP=="M12" & Chl.M$Chl_ug.cm2>6))),]

Chl.C.o<-Chl.C[-c(which((Chl.C$TimeP=="W2" & Chl.C$Chl_ug.cm2>2.5)| 
(Chl.C$TimeP=="M8" & Chl.C$Chl_ug.cm2>2.5) |
(Chl.C$TimeP=="M8" & Chl.C$Chl_ug.cm2<0.5))),]

Chl.H.o<-Chl.H[-c(which((Chl.H$TimeP=="W2" & Chl.H$Chl_ug.cm2>0.75) |
(Chl.H$TimeP=="M8" & Chl.H$Chl_ug.cm2>1.5) | 
(Chl.H$TimeP=="M12" & Chl.H$Chl_ug.cm2>0.51))),]

##Recombine Cleaned Chlorophyll dataframes
Chl.o<-rbind(Chl.M.o, Chl.C.o, Chl.H.o)

####Average Across Replicate Readings ("Rep" column: A, B, C)
names(Chl.o)
Chl.a<-aggregate(Chl.o$Chl_ug.cm2, list(Chl.o$RandN, Chl.o$ID), mean)
names(Chl.a)<-c("RandN", "ID", "Chl_ug.cm2")

####Add Total Chlorophyll to Main Coral Data and Thermal Tolerance Data

##Merge Averaged Chlorophyll Data with Coral Main Data
#Merges by Random Number (RandN) and ID columns
#Retains all Main samples
names(Chl.a)
CoralData<-merge(CoralData, Chl.a, all.x=TRUE, all.y=FALSE)

##Merge Averaged Chlorophyll Data with Thermal Tolerance Data
#Merges by Random Number (RandN) and ID columns
#Retains all Thermal Tolerance Assay samples
#Retains RandN, ID, TimeP, Site, Genotype, Treat, Treatment, Orig, Origin, Set, and SA_cm2 from Corals_Therm
names(Corals_Therm)
ThermData<-merge(Corals_Therm[,c(1:9, 13, 15)], Chl.a, all.x=TRUE, all.y=FALSE)


####Symbionts####

str(Sym)

####Calculate Dilution Factors

#Scales for 1:10 Dilution Sent for Flow Cytometry
Sym$Send_df<-Sym$SendInputVol_ul/Sym$SendTotalVol_ul

##Scales for Resuspending Symbiont Pellet based on Volume of 0.01% SDS Used
Sym$Resp_df<-Sym$InputVol_ul/Sym$RespVol_ul

####Calculate Symbiont Density

##Scale Flow Cy Count/ul to account for Dilution Factor applied at Flow Cy Center
Sym$Count_ul_scaled_FC_df<-Sym$Count_ul/Sym$FC_df

##Scale Flow Cy Count/ul to account for Sent Dilution Factor
Sym$Count_ul_scaled_send_df<-Sym$Count_ul_scaled_FC_df/Sym$Send_df

##Scale Flow Cy Count/ul to account for Resuspension Dilution Factor
Sym$Count_ul_scaled_resp_df<-Sym$Count_ul_scaled_send_df/Sym$Resp_df

##Merge with Sample Data to Calculate Symbionts per Surface Area
#Merges by Random Number (RandN) column
#Adds necessary Slurry Volume (Vol_ml) and Surface Area (SA_cm2) columns
#Retains Sample Meta Data only for samples with Symbiont data (Thermal Assay samples only)
Sym<-merge(Sym, SampData, all.x=TRUE, all.y=FALSE)

##Calculate Total Symbiont Count (scaling to Slurry Total Volume in ul)
Sym$Count_total<-Sym$Count_ul_scaled_resp_df*(Sym$Vol_ml*1000)

##Calculate Symbionts per Surface Area
Sym$Sym_cm2<-Sym$Count_total/Sym$SA_cm2

##Scale to Symb * 10^6 / cm^2
Sym$Sym10.6_cm2<-Sym$Sym_cm2/10^6

##Initial Visual Check
ggplot(Sym, aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))


####Check for Outliers

##Check for Outliers in Symbionts by Treatment and Timepoint

##Control
Sym.C<-subset(Sym, Treat=="C")
Sym.C<-Sym.C[-c(which(is.na(Sym.C$Sym10.6_cm2))),]

#Week 2
ggplot(subset(Sym.C, TimeP=="W2"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 90))

#Month 1
ggplot(subset(Sym.C, TimeP=="M1"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 90))

Sym.C$RandN[which(Sym.C$TimeP=="M1" & Sym.C$Sym10.6_cm2>1 & Sym.C$Site=="KL")] # "M1_14"

#Month 4
ggplot(subset(Sym.C, TimeP=="M4"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 90))

Sym.C$RandN[which(Sym.C$TimeP=="M4" & Sym.C$Sym10.6_cm2>0.85)] # "M4_53" "M4_53" "M4_74" "M4_76"
Sym.C$RandN[which(Sym.C$TimeP=="M4" & Sym.C$Sym10.6_cm2>0.45 & Sym.C$Site=="KL")] # "M4_29" "M4_29" "M4_83" "M4_83"

#Month 8
ggplot(subset(Sym.C, TimeP=="M8"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 90))

Sym.C$RandN[which(Sym.C$TimeP=="M8" & Sym.C$Sym10.6_cm2>1.5)] #"M8_33" "M8_33" "M8_33"
Sym.C$RandN[which(Sym.C$TimeP=="M8" & Sym.C$Sym10.6_cm2>1 & Sym.C$Site=="KL")] #"M8_48"

#Month 12
ggplot(subset(Sym.C, TimeP=="M12"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 90))


##Heated
Sym.H<-subset(Sym, Treat=="H")

#Week 2
ggplot(subset(Sym.H, TimeP=="W2"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90))

#Month 1
ggplot(subset(Sym.H, TimeP=="M1"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90))

Sym.H$RandN[which(Sym.H$TimeP=="M1" & Sym.H$Sym10.6_cm2>0.75)] #"M1_73" "M1_73"
Sym.H$RandN[which(Sym.H$TimeP=="M1" & Sym.H$Sym10.6_cm2>0.55 & Sym.H$Site=="KL")] #"M1_58" "M1_58"
Sym.H$RandN[which(Sym.H$TimeP=="M1" & Sym.H$Sym10.6_cm2>0.4 & Sym.H$Site=="SS" & Sym.H$Genotype=="AC8")] # "M1_94" "M1_94" "M1_94"

#Month 4
ggplot(subset(Sym.H, TimeP=="M4"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1.5)+
  theme(axis.text.x = element_text(angle = 90))

Sym.H$RandN[which(Sym.H$TimeP=="M4" & Sym.H$Sym10.6_cm2>0.6)] #"M4_44" "M4_49" "M4_61" "M4_61" "M4_61" "M4_85"
Sym.H$RandN[which(Sym.H$TimeP=="M4" & Sym.H$Sym10.6_cm2>0.45 & Sym.H$Site=="KL")] #"M4_59" "M4_72" "M4_72" "M4_88" "M4_88" "M4_94"

#Month 8
ggplot(subset(Sym.H, TimeP=="M8"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,2.1)+
  theme(axis.text.x = element_text(angle = 90))

Sym.H$RandN[which(Sym.H$TimeP=="M8" & Sym.H$Sym10.6_cm2>1.5)] #"M8_42" "M8_42" "M8_42"
Sym.H$RandN[which(Sym.H$TimeP=="M8" & Sym.H$Sym10.6_cm2>1.2 & Sym.H$Site=="KL")] #"M8_28" "M8_28" "M8_28"

#Month 12
ggplot(subset(Sym.H, TimeP=="M12"), aes(x=Set, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 90))

Sym.H$RandN[which(Sym.H$TimeP=="M12" & Sym.H$Sym10.6_cm2>0.5)] #"M12_48" "M12_48" "M12_48" "M12_63" "M12_63" "M12_63"



####Fv/Fm####
str(PAM)
str(PAM_Meta)

####Organize Data

##Merge PAM Meta Data with PAM Data
#Merges by Memory Number (Memory) column
#Retains PAM Data only for relevant measurements with Meta Data 
PAM<-merge(PAM_Meta, PAM, all.x=TRUE, all.y=FALSE)

PAM[(which(is.na(PAM$Fv_Fm))),]

##Merge PAM Data with Sample Data
#Merges by ID column
#Retains Sample Meta Data only for samples with Fv/Fm data (Thermal Assay samples only)
PAM<-merge(PAM, SampData, all.x=TRUE, all.y=FALSE)

##Initial Visual Check
ggplot(PAM, aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))


####Check for Outliers

##Check for Outliers in Symbionts by Treatment and Timepoint

##Control
PAM.C<-subset(PAM, Treat=="C")
PAM.C<-PAM.C[-c(which(is.na(PAM.C$Fv_Fm))),]

#Week 2
ggplot(subset(PAM.C, TimeP=="W2"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.45,0.75)+
  theme(axis.text.x = element_text(angle = 90))

PAM.C$RandN[which(PAM.C$TimeP=="W2" & PAM.C$Fv_Fm<0.5)] # "W2_93" "W2_71"

#Month 1
ggplot(subset(PAM.C, TimeP=="M1"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.45,0.75)+
  theme(axis.text.x = element_text(angle = 90))

#Month 4
ggplot(subset(PAM.C, TimeP=="M4"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.45,0.75)+
  theme(axis.text.x = element_text(angle = 90))

PAM.C$RandN[which(PAM.C$TimeP=="M4" & PAM.C$Fv_Fm<0.53)] # "M4_91"

#Month 8
ggplot(subset(PAM.C, TimeP=="M8"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.45,0.75)+
  theme(axis.text.x = element_text(angle = 90))

#Month 12
ggplot(subset(PAM.C, TimeP=="M12"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.45,0.75)+
  theme(axis.text.x = element_text(angle = 90))


##Heated
PAM.H<-subset(PAM, Treat=="H")
PAM.H<-PAM.H[-c(which(is.na(PAM.H$Fv_Fm))),]

#Week 2
ggplot(subset(PAM.H, TimeP=="W2"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.25,0.7)+
  theme(axis.text.x = element_text(angle = 90))

PAM.H$RandN[which(PAM.H$TimeP=="W2" & PAM.H$Fv_Fm<0.5 & PAM.H$Site=="KL" & PAM.H$Genotype=="AC10")] # "W2_89"
PAM.H$RandN[which(PAM.H$TimeP=="W2" & PAM.H$Fv_Fm<0.5 & PAM.H$Site=="SS" & PAM.H$Genotype=="AC12" & PAM.H$Orig=="N")] # "W2_94" "W2_94"
PAM.H$RandN[which(PAM.H$TimeP=="W2" & PAM.H$Fv_Fm<0.5 & PAM.H$Site=="KL" & PAM.H$Genotype=="AC12" & PAM.H$Orig=="T")] # "W2_38" "W2_38" "W2_38"

#Month 1
ggplot(subset(PAM.H, TimeP=="M1"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.25,0.7)+
  theme(axis.text.x = element_text(angle = 90))

PAM.H$RandN[which(PAM.H$TimeP=="M1" & PAM.H$Fv_Fm<0.35)] # "M1_86"

#Month 4
ggplot(subset(PAM.H, TimeP=="M4"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.25,0.7)+
  theme(axis.text.x = element_text(angle = 90))

PAM.H$RandN[which(PAM.H$TimeP=="M4" & PAM.H$Fv_Fm<0.35 &  PAM.H$Site=="SS")] # "M4_84"

#Month 8
ggplot(subset(PAM.H, TimeP=="M8"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.25,0.7)+
  theme(axis.text.x = element_text(angle = 90))

PAM.H$RandN[which(PAM.H$TimeP=="M8" & PAM.H$Fv_Fm<0.45)] # "M8_83"

#Month 12
ggplot(subset(PAM.H, TimeP=="M12"), aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  ylim(0.2,0.7)+
  theme(axis.text.x = element_text(angle = 90))


##Remove Outliers 
PAM.C.o<-PAM.C[-c(which((PAM.C$TimeP=="W2" & PAM.C$Fv_Fm<0.5) |
              (PAM.C$TimeP=="M4" & PAM.C$Fv_Fm<0.53))),]

PAM.H.o<-PAM.H[-c(which((PAM.H$TimeP=="W2" & PAM.H$Fv_Fm<0.5 & PAM.H$Site=="KL" & PAM.H$Genotype=="AC10") | 
              (PAM.H$TimeP=="W2" & PAM.H$Fv_Fm<0.5 & PAM.H$Site=="SS" & PAM.H$Genotype=="AC12" & PAM.H$Orig=="N") |
              (PAM.H$TimeP=="W2" & PAM.H$Fv_Fm<0.5 & PAM.H$Site=="KL" & PAM.H$Genotype=="AC12" & PAM.H$Orig=="T") |
              (PAM.H$TimeP=="M1" & PAM.H$Fv_Fm<0.35)|
              (PAM.H$TimeP=="M4" & PAM.H$Fv_Fm<0.35 &  PAM.H$Site=="SS") |
                (PAM.H$TimeP=="M8" & PAM.H$Fv_Fm<0.45))),]

##Recombine Cleaned PAM dataframes
PAM.o<-rbind(PAM.C.o, PAM.H.o)

ggplot(PAM.o, aes(x=Set, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(axis.text.x = element_text(angle = 90))


####Average Across Replicate Readings ("Rep" column: A, B, C)
names(PAM.o)
PAM.a<-aggregate(PAM.o$Fv_Fm, list(PAM.o$ID), mean)
names(PAM.a)<-c("ID", "Fv_Fm")

####Add Fv/Fm to Thermal Tolerance Data

##Merge Averaged Fv/Fm Data with Thermal Tolerance Data
#Merges by ID column
#Retains all Thermal Tolerance Assay samples
names(ThermData)
ThermData<-merge(ThermData, PAM.a, all.x=TRUE, all.y=FALSE)
