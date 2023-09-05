#Title: "Coral Thermal Tolerance following Reciprocal Transplant"
#Author: "Serena Hackerott"
#Date: "08/29/2023"

#-------Set Up------------

####Load Packages####
library(Rmisc)
library(ggplot2)

#Note: Run "Graphing Parameters" section from 01_ExperimentalSetup.R file


####Load Data####
#Note: Physiological metrics calculated in 02_PhysiologyMetrics.R file
Thermal<-read.csv("Outputs/ThermData.csv", header=TRUE)


#-------Confirm Bleaching------------

####Week 2####
Therm.W2<-subset(Thermal, TimeP=="W2")

####Chlorophyll

##Check normality
hist(Therm.W2$Chl_ug.cm2)
qqnorm(Therm.W2$Chl_ug.cm2)
shapiro.test(Therm.W2$Chl_ug.cm2)
#Not Normal

hist(log(Therm.W2$Chl_ug.cm2+1))
qqnorm(log(Therm.W2$Chl_ug.cm2+1))
shapiro.test(log(Therm.W2$Chl_ug.cm2+1))
##Not Normal

kruskal.test(Chl_ug.cm2~Treatment, Therm.W2)
# Kruskal-Wallis chi-squared = 53.055, df = 1, p-value = 3.244e-13

ggplot(Therm.W2, aes(x=Treatment, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Symbionts

##Check normality
hist(Therm.W2$Sym10.6_cm2)
qqnorm(Therm.W2$Sym10.6_cm2)
shapiro.test(Therm.W2$Sym10.6_cm2)
#Not Normal

hist(log(Therm.W2$Sym10.6_cm2+1))
qqnorm(log(Therm.W2$Sym10.6_cm2+1))
shapiro.test(log(Therm.W2$Sym10.6_cm2+1))
##Not Normal

kruskal.test(Sym10.6_cm2~Treatment, Therm.W2)
# Kruskal-Wallis chi-squared = 27.681, df = 1, p-value = 1.431e-07

ggplot(Therm.W2, aes(x=Treatment, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Fv/Fm

##Check normality
hist(Therm.W2$Fv_Fm)
qqnorm(Therm.W2$Fv_Fm)
shapiro.test(Therm.W2$Fv_Fm)
#Not Normal

hist(log(Therm.W2$Fv_Fm+1))
qqnorm(log(Therm.W2$Fv_Fm+1))
shapiro.test(log(Therm.W2$Fv_Fm+1))
##Not Normal

kruskal.test(Fv_Fm~Treatment, Therm.W2)
# Kruskal-Wallis chi-squared = 46.444, df = 1, p-value = 9.425e-12

ggplot(Therm.W2, aes(x=Treatment, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))



####Month 1####
Therm.M1<-subset(Thermal, TimeP=="M1")

####Chlorophyll

##Check normality
hist(Therm.M1$Chl_ug.cm2)
qqnorm(Therm.M1$Chl_ug.cm2)
shapiro.test(Therm.M1$Chl_ug.cm2)
#Not Normal

hist(log(Therm.M1$Chl_ug.cm2+1))
qqnorm(log(Therm.M1$Chl_ug.cm2+1))
shapiro.test(log(Therm.M1$Chl_ug.cm2+1))
##Not Normal

kruskal.test(Chl_ug.cm2~Treatment, Therm.M1)
# Kruskal-Wallis chi-squared = 71.258, df = 1, p-value < 2.2e-16

ggplot(Therm.M1, aes(x=Treatment, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Symbionts

##Check normality
hist(Therm.M1$Sym10.6_cm2)
qqnorm(Therm.M1$Sym10.6_cm2)
shapiro.test(Therm.M1$Sym10.6_cm2)
#Not Normal

hist(log(Therm.M1$Sym10.6_cm2+1))
qqnorm(log(Therm.M1$Sym10.6_cm2+1))
shapiro.test(log(Therm.M1$Sym10.6_cm2+1))
##Not Normal

kruskal.test(Sym10.6_cm2~Treatment, Therm.M1)
# Kruskal-Wallis chi-squared = 41.77, df = 1, p-value = 1.027e-10

ggplot(Therm.M1, aes(x=Treatment, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Fv/Fm

##Check normality
hist(Therm.M1$Fv_Fm)
qqnorm(Therm.M1$Fv_Fm)
shapiro.test(Therm.M1$Fv_Fm)
#Not Normal

hist(log(Therm.M1$Fv_Fm+1))
qqnorm(log(Therm.M1$Fv_Fm+1))
shapiro.test(log(Therm.M1$Fv_Fm+1))
##Not Normal

kruskal.test(Fv_Fm~Treatment, Therm.M1)
# Kruskal-Wallis chi-squared = 62.863, df = 1, p-value = 2.216e-15

ggplot(Therm.M1, aes(x=Treatment, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Month 4####
Therm.M4<-subset(Thermal, TimeP=="M4")

####Chlorophyll

##Check normality
hist(Therm.M4$Chl_ug.cm2)
qqnorm(Therm.M4$Chl_ug.cm2)
shapiro.test(Therm.M4$Chl_ug.cm2)
#Not Normal

hist(log(Therm.M4$Chl_ug.cm2+1))
qqnorm(log(Therm.M4$Chl_ug.cm2+1))
shapiro.test(log(Therm.M4$Chl_ug.cm2+1))
##Not Normal

kruskal.test(Chl_ug.cm2~Treatment, Therm.M4)
# Kruskal-Wallis chi-squared = 64.97, df = 1, p-value = 7.605e-16

ggplot(Therm.M4, aes(x=Treatment, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Symbionts

##Check normality
hist(Therm.M4$Sym10.6_cm2)
qqnorm(Therm.M4$Sym10.6_cm2)
shapiro.test(Therm.M4$Sym10.6_cm2)
#Not Normal

hist(log(Therm.M4$Sym10.6_cm2+1))
qqnorm(log(Therm.M4$Sym10.6_cm2+1))
shapiro.test(log(Therm.M4$Sym10.6_cm2+1))
##Not Normal

kruskal.test(Sym10.6_cm2~Treatment, Therm.M4)
# Kruskal-Wallis chi-squared = 6.5026, df = 1, p-value = 0.01077

ggplot(Therm.M4, aes(x=Treatment, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Fv/Fm

##Check normality
hist(Therm.M4$Fv_Fm)
qqnorm(Therm.M4$Fv_Fm)
shapiro.test(Therm.M4$Fv_Fm)
#Not Normal

hist(log(Therm.M4$Fv_Fm+1))
qqnorm(log(Therm.M4$Fv_Fm+1))
shapiro.test(log(Therm.M4$Fv_Fm+1))
##Not Normal

kruskal.test(Fv_Fm~Treatment, Therm.M4)
# Kruskal-Wallis chi-squared = 59.655, df = 1, p-value = 1.13e-14

ggplot(Therm.M4, aes(x=Treatment, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))

####Month 8####
Therm.M8<-subset(Thermal, TimeP=="M8")

####Chlorophyll

##Check normality
hist(Therm.M8$Chl_ug.cm2)
qqnorm(Therm.M8$Chl_ug.cm2)
shapiro.test(Therm.M8$Chl_ug.cm2)
#Not Normal

hist(log(Therm.M8$Chl_ug.cm2+1))
qqnorm(log(Therm.M8$Chl_ug.cm2+1))
shapiro.test(log(Therm.M8$Chl_ug.cm2+1))
##Not Normal

kruskal.test(Chl_ug.cm2~Treatment, Therm.M8)
# Kruskal-Wallis chi-squared = 25.696, df = 1, p-value = 3.997e-07

ggplot(Therm.M8, aes(x=Treatment, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Symbionts

##Check normality
hist(Therm.M8$Sym10.6_cm2)
qqnorm(Therm.M8$Sym10.6_cm2)
shapiro.test(Therm.M8$Sym10.6_cm2)
#Not Normal

hist(log(Therm.M8$Sym10.6_cm2+1))
qqnorm(log(Therm.M8$Sym10.6_cm2+1))
shapiro.test(log(Therm.M8$Sym10.6_cm2+1))
##Not Normal

kruskal.test(Sym10.6_cm2~Treatment, Therm.M8)
# Kruskal-Wallis chi-squared = 1.0078, df = 1, p-value = 0.3154

ggplot(Therm.M8, aes(x=Treatment, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Fv/Fm

##Check normality
hist(Therm.M8$Fv_Fm)
qqnorm(Therm.M8$Fv_Fm)
shapiro.test(Therm.M8$Fv_Fm)
#Normal

bartlett.test(Fv_Fm~Treatment, Therm.M8)
#Equal Variance 

summary(aov(Fv_Fm~Treatment, Therm.M8))
# Df  Sum Sq  Mean Sq F value   Pr(>F)    
# Treatment    1 0.01529 0.015293   30.67 2.76e-07 ***
# Residuals   94 0.04687 0.000499 

ggplot(Therm.M8, aes(x=Treatment, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Month 12####
Therm.M12<-subset(Thermal, TimeP=="M12")

####Chlorophyll

##Check normality
hist(Therm.M12$Chl_ug.cm2)
qqnorm(Therm.M12$Chl_ug.cm2)
shapiro.test(Therm.M12$Chl_ug.cm2)
#Not Normal

hist(log(Therm.M12$Chl_ug.cm2+1))
qqnorm(log(Therm.M12$Chl_ug.cm2+1))
shapiro.test(log(Therm.M12$Chl_ug.cm2+1))
##Not Normal

kruskal.test(Chl_ug.cm2~Treatment, Therm.M12)
# Kruskal-Wallis chi-squared = 70.5, df = 1, p-value < 2.2e-16

ggplot(Therm.M12, aes(x=Treatment, y=Chl_ug.cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Symbionts

##Check normality
hist(Therm.M12$Sym10.6_cm2)
qqnorm(Therm.M12$Sym10.6_cm2)
shapiro.test(Therm.M12$Sym10.6_cm2)
#Not Normal

hist(log(Therm.M12$Sym10.6_cm2+1))
qqnorm(log(Therm.M12$Sym10.6_cm2+1))
shapiro.test(log(Therm.M12$Sym10.6_cm2+1))
##Not Normal

kruskal.test(Sym10.6_cm2~Treatment, Therm.M12)
# Kruskal-Wallis chi-squared = 65.325, df = 1, p-value = 6.352e-16

ggplot(Therm.M12, aes(x=Treatment, y=Sym10.6_cm2)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))


####Fv/Fm

##Check normality
hist(Therm.M12$Fv_Fm)
qqnorm(Therm.M12$Fv_Fm)
shapiro.test(Therm.M12$Fv_Fm)
#Not Normal

hist(log(Therm.M12$Fv_Fm+1))
qqnorm(log(Therm.M12$Fv_Fm+1))
shapiro.test(log(Therm.M12$Fv_Fm+1))
##Not Normal

kruskal.test(Fv_Fm~Treatment, Therm.M12)
# Kruskal-Wallis chi-squared = 71.261, df = 1, p-value < 2.2e-16

ggplot(Therm.M12, aes(x=Treatment, y=Fv_Fm)) + 
  geom_boxplot(alpha=0.5, shape=2, outlier.shape = NA)+
  geom_jitter(shape=16, position=position_jitter(0.1))



#-------Thermal Tolerance (Retention)------------

####Calculate Retention####

####Averages by Group for Controls

##Subset Controls
Therm.C<-subset(Thermal, Treat=="C")
names(Therm.C)

##Calculate Averages of Chlorophyll, Symbionts, and Fv/Fm by Timepoint, Genotype, Site, and Origin
Therm.C.a<-aggregate(Therm.C[,c(12:14)], list(Therm.C$TimeP, Therm.C$Site, Therm.C$Genotype, Therm.C$Orig), mean, na.action=na.omit)
names(Therm.C.a)[1:4]<-c("TimeP", "Site", "Genotype", "Orig")
names(Therm.C.a)[5:7]<-paste(names(Therm.C.a)[5:7], "C", sep="_")

##Fix one NA averaging issue
Therm.C.a$Chl_ug.cm2_C[which(Therm.C.a$TimeP=="M8" & Therm.C.a$Site=="SS" & 
                               Therm.C.a$Orig== "T" &Therm.C.a$Genotype=="AC12")]<-
  mean(na.omit(Therm.C$Chl_ug.cm2[which(Therm.C$TimeP=="M8" & Therm.C$Site=="SS" & 
  Therm.C$Orig== "T" & Therm.C$Genotype=="AC12")]))


####Calculate Retention of Heated Relative to Controls

##Subset Heated
Therm.H<-subset(Thermal, Treat=="H")
names(Therm.H)

##Merge Heated Dataframe with Control Averages
#Retains all Heated samples
#Merges by Timepoint, Genotype, Site, and Origin
#Adds average Control Chlorophyll, Symbionts, and Fv/Fm values
Therm.H<-merge(Therm.H, Therm.C.a, all.x=TRUE, all.y=FALSE)

##Calculate Percent Retention of each Bleaching Metric
Therm.H$Chl_p<-(Therm.H$Chl_ug.cm2/Therm.H$Chl_ug.cm2_C)*100
Therm.H$Sym_p<-(Therm.H$Sym10.6_cm2/Therm.H$Sym10.6_cm2_C)*100
Therm.H$Fv_Fm_p<-(Therm.H$Fv_Fm/Therm.H$Fv_Fm_C)*100

##Set Retention >100% to 100%
Therm.H$Chl_p[which(Therm.H$Chl_p>100)]<-100
Therm.H$Sym_p[which(Therm.H$Sym_p>100)]<-100
Therm.H$Fv_Fm_p[which(Therm.H$Fv_Fm_p>100)]<-100



####Week 2####
Therm.H.W2<-subset(Therm.H, TimeP=="W2")

####Chlorophyll

##Check normality
hist(Therm.H.W2$Chl_p)
qqnorm(Therm.H.W2$Chl_p)
shapiro.test(Therm.H.W2$Chl_p)
#Normal


####Symbionts

##Check normality
hist(Therm.H.W2$Sym_p)
qqnorm(Therm.H.W2$Sym_p)
shapiro.test(Therm.H.W2$Sym_p)
##Not Normal


####Fv/Fm

##Check normality
hist(Therm.H.W2$Fv_Fm_p)
qqnorm(Therm.H.W2$Fv_Fm_p)
shapiro.test(Therm.H.W2$Fv_Fm_p)
##Not Normal


####Figures####

####Summary Statistics 

##Chlorophyll

Chl_p.SE<-summarySE(Therm.H, measurevar="Chl_p", groupvars=c("TimeP", "Genotype", "Site", "Orig", "Origin"), na.rm=TRUE)

str(Chl_p.SE)
Chl_p.SE$TimeP<-factor(Chl_p.SE$TimeP, levels=c("W2", "M1", "M4", "M8", "M12"), ordered=TRUE)
Chl_p.SE$Site<-factor(Chl_p.SE$Site, levels=c("SS", "KL"), ordered=TRUE)
Chl_p.SE$Genotype<-factor(Chl_p.SE$Genotype, levels=c("AC8", "AC10", "AC12"), ordered=TRUE)
Chl_p.SE$Orig<-factor(Chl_p.SE$Orig, levels=c("N", "T"), ordered=TRUE)
Chl_p.SE$Origin<-factor(Chl_p.SE$Origin, levels=c("Native", "Transplant"), ordered=TRUE)

Chl_p.SE$Geno.Orig<-paste(Chl_p.SE$Genotype, Chl_p.SE$Orig, sep=".")
Chl_p.SE$Geno.Site<-paste(Chl_p.SE$Genotype, Chl_p.SE$Site,  sep=".")
Chl_p.SE$Geno.Site<-factor(Chl_p.SE$Geno.Site, levels=c("AC8.SS", "AC8.KL", "AC10.SS",  "AC10.KL", "AC12.SS",  "AC12.KL"), ordered=TRUE)



####Plot Mean +/- Standard Deviation by Genotype and Site

##Week 2
Chl_p.SE.W2.plot<-ggplot(Chl_p.SE[which(Chl_p.SE$TimeP=="W2"),], aes(x=Geno.Site, y=Chl_p, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_p-sd, ymax=Chl_p+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 50)+
  labs(x="Genotype and Site", y=expression(paste('% Retention of Chlorophyll (\u03BCg cm'^-2*")")));Chl_p.SE.W2.plot


##Month 1
Chl_p.SE.M1.plot<-ggplot(Chl_p.SE[which(Chl_p.SE$TimeP=="M1"),], aes(x=Geno.Site, y=Chl_p, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_p-sd, ymax=Chl_p+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 50)+
  labs(x="Genotype and Site", y=expression(paste('% Retention of Chlorophyll (\u03BCg cm'^-2*")")));Chl_p.SE.M1.plot


##Month 4
Chl_p.SE.M4.plot<-ggplot(Chl_p.SE[which(Chl_p.SE$TimeP=="M4"),], aes(x=Geno.Site, y=Chl_p, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_p-sd, ymax=Chl_p+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 50)+
  labs(x="Genotype and Site", y=expression(paste('% Retention of Chlorophyll (\u03BCg cm'^-2*")")));Chl_p.SE.M4.plot


##Month 8
Chl_p.SE.M8.plot<-ggplot(Chl_p.SE[which(Chl_p.SE$TimeP=="M8"),], aes(x=Geno.Site, y=Chl_p, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_p-sd, ymax=Chl_p+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 110)+
  labs(x="Genotype and Site", y=expression(paste('% Retention of Chlorophyll (\u03BCg cm'^-2*")")));Chl_p.SE.M8.plot


##Month 12
Chl_p.SE.M12.plot<-ggplot(Chl_p.SE[which(Chl_p.SE$TimeP=="M12"),], aes(x=Geno.Site, y=Chl_p, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_p-sd, ymax=Chl_p+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 50)+
  labs(x="Genotype and Site", y=expression(paste('% Retention of Chlorophyll (\u03BCg cm'^-2*")")));Chl_p.SE.M12.plot
