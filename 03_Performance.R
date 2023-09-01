#Title: "Coral Performance following Reciprocal Transplant"
#Author: "Serena Hackerott"
#Date: "08/29/2023"

#-------Set Up------------

####Load Packages####
library(ggplot2)
library(effectsize)
library(emmeans)

#Note: Run "Graphing Parameters" section from 01_ExperimentalSetup.R file


####Load Data####
#Note: Physiological metrics calculated in 02_PhysiologyMetrics.R file
Coral<-read.csv("Outputs/CoralData.csv", header=TRUE)


#-------Univariate Physiology------------

####Protein Host Fraction####

##Check normality
hist(Coral$TP_ug.cm2_C)
qqnorm(Coral$TP_ug.cm2_C)
shapiro.test(Coral$TP_ug.cm2_C)
#Not normal

hist(log(Coral$TP_ug.cm2_C+1))
qqnorm(log(Coral$TP_ug.cm2_C+1))
shapiro.test(log(Coral$TP_ug.cm2_C+1))
#Not Normal but less Skewed

##Model with log+1 transformation
Prot.C.lm<-lm(log(TP_ug.cm2_C+1)~Genotype+Site+Origin+TimeP, data=Coral)

##Check residuals
plot(fitted(Prot.C.lm), resid(Prot.C.lm))
abline(0,0)

qqnorm(resid(Prot.C.lm))
qqline(resid(Prot.C.lm))

plot(density(resid(Prot.C.lm)))

##Model results
summary(Prot.C.lm)
# Multiple R-squared:  0.2618,	Adjusted R-squared:  0.2359 
# F-statistic: 10.11 on 8 and 228 DF,  p-value: 4.699e-12

anova(Prot.C.lm)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype    2 1.2434 0.62168 21.3443 3.184e-09 ***
# Site        1 0.7154 0.71543 24.5629 1.403e-06 ***
# Origin      1 0.0334 0.03337  1.1458   0.28557    
# TimeP       4 0.3632 0.09081  3.1176   0.01596 *  
# Residuals 228 6.6408 0.02913  

eta_squared(Prot.C.lm, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype  |     0.14 | [0.07, 1.00]
# Site      |     0.08 | [0.03, 1.00]
# Origin    | 3.71e-03 | [0.00, 1.00]
# TimeP     |     0.04 | [0.00, 1.00]

emmeans(Prot.C.lm, pairwise~TimeP)
# contrast  estimate     SE  df t.ratio p.value
# M1 - M12  9.22e-02 0.0354 228   2.602  0.0733
# M1 - M4   1.04e-01 0.0350 228   2.960  0.0278 *
# M1 - M8   9.21e-02 0.0350 228   2.630  0.0683
# M1 - W2   4.26e-02 0.0350 228   1.216  0.7421
# M12 - M4  1.15e-02 0.0352 228   0.327  0.9975
# M12 - M8 -4.77e-05 0.0352 228  -0.001  1.0000
# M12 - W2 -4.96e-02 0.0352 228  -1.407  0.6236
# M4 - M8  -1.16e-02 0.0348 228  -0.332  0.9974
# M4 - W2  -6.11e-02 0.0348 228  -1.754  0.4034
# M8 - W2  -4.95e-02 0.0348 228  -1.422  0.6144


####Protein Symbiont Fraction####

##Check normality
hist(Coral$TP_ug.cm2_S)
qqnorm(Coral$TP_ug.cm2_S)
shapiro.test(Coral$TP_ug.cm2_S)
#Not normal

hist(log(Coral$TP_ug.cm2_S+1))
qqnorm(log(Coral$TP_ug.cm2_S+1))
shapiro.test(log(Coral$TP_ug.cm2_S+1))
#Normal

##Model with log+1 transformation
Prot.S.lm<-lm(log(TP_ug.cm2_S+1)~Genotype+Site+Origin+TimeP, data=Coral)

##Check residuals
plot(fitted(Prot.S.lm), resid(Prot.S.lm))
abline(0,0)

qqnorm(resid(Prot.S.lm))
qqline(resid(Prot.S.lm))

plot(density(resid(Prot.S.lm)))

##Model results
summary(Prot.S.lm)
# Multiple R-squared:  0.5555,	Adjusted R-squared:   0.54 
# F-statistic: 35.78 on 8 and 229 DF,  p-value: < 2.2e-16

anova(Prot.S.lm)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype    2 0.9439 0.47197  15.793 3.758e-07 ***
# Site        1 1.0434 1.04341  34.913 1.239e-08 ***
# Origin      1 0.0633 0.06327   2.117     0.147    
# TimeP       4 6.5032 1.62579  54.400 < 2.2e-16 ***
# Residuals 229 6.8438 0.02989   

eta_squared(Prot.S.lm, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype  |     0.06 | [0.02, 1.00]
# Site      |     0.07 | [0.02, 1.00]
# Origin    | 4.11e-03 | [0.00, 1.00]
# TimeP     |     0.42 | [0.34, 1.00]

emmeans(Prot.S.lm, pairwise~TimeP)
# contrast estimate     SE  df t.ratio p.value
# M1 - M12   0.4270 0.0357 229  11.972  <.0001 *
# M1 - M4    0.3606 0.0355 229  10.165  <.0001 *
# M1 - M8    0.2855 0.0355 229   8.047  <.0001 *
# M1 - W2    0.0701 0.0355 229   1.975  0.2813
# M12 - M4  -0.0664 0.0355 229  -1.871  0.3360
# M12 - M8  -0.1415 0.0355 229  -3.989  0.0008 *
# M12 - W2  -0.3569 0.0355 229 -10.061  <.0001 *
# M4 - M8   -0.0751 0.0353 229  -2.129  0.2113
# M4 - W2   -0.2906 0.0353 229  -8.234  <.0001 *
# M8 - W2   -0.2154 0.0353 229  -6.105  <.0001 *


####Biomass Host Fraction####

##Check normality
hist(Coral$AFDW_mg.cm2_C)
qqnorm(Coral$AFDW_mg.cm2_C)
shapiro.test(Coral$AFDW_mg.cm2_C)
#Not normal

hist(log(Coral$AFDW_mg.cm2_C+1))
qqnorm(log(Coral$AFDW_mg.cm2_C+1))
shapiro.test(log(Coral$AFDW_mg.cm2_C+1))
#Not Normal but less Skewed

##Model with log+1 transformation
Bio.C.lm<-lm(log(AFDW_mg.cm2_C+1)~Genotype+Site+Origin+TimeP, data=Coral)

##Check residuals
plot(fitted(Bio.C.lm), resid(Bio.C.lm))
abline(0,0)

qqnorm(resid(Bio.C.lm))
qqline(resid(Bio.C.lm))

plot(density(resid(Bio.C.lm)))

##Model results
summary(Bio.C.lm)
# Multiple R-squared:  0.1529,	Adjusted R-squared:  0.1233 
# F-statistic: 5.166 on 8 and 229 DF,  p-value: 6.12e-06

anova(Bio.C.lm)
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Genotype    2 0.2860 0.143002 10.1524 5.963e-05 ***
# Site        1 0.1085 0.108503  7.7032  0.005968 ** 
# Origin      1 0.0179 0.017867  1.2684  0.261238    
# TimeP       4 0.1698 0.042449  3.0137  0.018909 *  
# Residuals 229 3.2256 0.014085     

eta_squared(Bio.C.lm, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype  |     0.08 | [0.03, 1.00]
# Site      |     0.03 | [0.00, 1.00]
# Origin    | 4.69e-03 | [0.00, 1.00]
# TimeP     |     0.04 | [0.00, 1.00]

emmeans(Bio.C.lm, pairwise~TimeP)
# contrast estimate     SE  df t.ratio p.value
# M1 - M12  0.06005 0.0245 229   2.452  0.1054
# M1 - M4   0.06319 0.0244 229   2.594  0.0748
# M1 - M8   0.02690 0.0244 229   1.105  0.8040
# M1 - W2   0.07103 0.0244 229   2.916  0.0316 *
# M12 - M4  0.00313 0.0244 229   0.129  0.9999
# M12 - M8 -0.03315 0.0244 229  -1.361  0.6532
# M12 - W2  0.01097 0.0244 229   0.451  0.9914
# M4 - M8  -0.03628 0.0242 229  -1.498  0.5651
# M4 - W2   0.00784 0.0242 229   0.324  0.9976
# M8 - W2   0.04412 0.0242 229   1.821  0.3638


####Biomass Symbiont Fraction####

##Check normality
hist(Coral$AFDW_mg.cm2_S)
qqnorm(Coral$AFDW_mg.cm2_S)
shapiro.test(Coral$AFDW_mg.cm2_S)
#Normal

##Model 
Bio.S.lm<-lm(AFDW_mg.cm2_S~Genotype+Site+Origin+TimeP, data=Coral)

##Check residuals
plot(fitted(Bio.S.lm), resid(Bio.S.lm))
abline(0,0)

qqnorm(resid(Bio.S.lm))
qqline(resid(Bio.S.lm))

plot(density(resid(Bio.S.lm)))

##Model results
summary(Bio.S.lm)
# Multiple R-squared:  0.2182,	Adjusted R-squared:  0.1909 
# F-statistic: 7.988 on 8 and 229 DF,  p-value: 1.733e-09

anova(Bio.S.lm)
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Genotype    2 0.2841 0.14204  4.7683 0.0093566 ** 
# Site        1 0.8900 0.89001 29.8780 1.199e-07 ***
# Origin      1 0.0438 0.04379  1.4702 0.2265699    
# TimeP       4 0.6858 0.17145  5.7557 0.0001974 ***
# Residuals 229 6.8215 0.02979        

eta_squared(Bio.S.lm, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype  |     0.03 | [0.00, 1.00]
# Site      |     0.10 | [0.05, 1.00]
# Origin    | 5.02e-03 | [0.00, 1.00]
# TimeP     |     0.08 | [0.02, 1.00]

emmeans(Bio.S.lm, pairwise~TimeP)
# contrast estimate     SE  df t.ratio p.value
# M1 - M12  0.06416 0.0356 229   1.802  0.3751
# M1 - M4  -0.02682 0.0354 229  -0.757  0.9424
# M1 - M8  -0.00968 0.0354 229  -0.273  0.9988
# M1 - W2   0.11642 0.0354 229   3.287  0.0102 *
# M12 - M4 -0.09098 0.0354 229  -2.569  0.0797
# M12 - M8 -0.07383 0.0354 229  -2.084  0.2303
# M12 - W2  0.05227 0.0354 229   1.476  0.5794
# M4 - M8   0.01714 0.0352 229   0.487  0.9885
# M4 - W2   0.14324 0.0352 229   4.066  0.0006 *
# M8 - W2   0.12610 0.0352 229   3.579  0.0038


####Chlorophyll####

##Check normality
hist(Coral$Chl_ug.cm2)
qqnorm(Coral$Chl_ug.cm2)
shapiro.test(Coral$Chl_ug.cm2)
#Not normal

hist(log(Coral$Chl_ug.cm2+1))
qqnorm(log(Coral$Chl_ug.cm2+1))
shapiro.test(log(Coral$Chl_ug.cm2+1))
#Not Normal but less Skewed

##Model with log+1 transformation
Chl.lm<-lm(log(Chl_ug.cm2+1)~Genotype+Site+Origin+TimeP, data=Coral)

##Check residuals
plot(fitted(Chl.lm), resid(Chl.lm))
abline(0,0)

qqnorm(resid(Chl.lm))
qqline(resid(Chl.lm))

plot(density(resid(Chl.lm)))

##Model results
summary(Chl.lm)
# Multiple R-squared:  0.6802,	Adjusted R-squared:  0.669 
# F-statistic: 60.88 on 8 and 229 DF,  p-value: < 2.2e-16

anova(Chl.lm)
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Genotype    2 2.6400  1.3200  62.7816 <2e-16 ***
# Site        1 3.5074  3.5074 166.8174 <2e-16 ***
# Origin      1 0.0480  0.0480   2.2852  0.132    
# TimeP       4 4.0453  1.0113  48.1008 <2e-16 ***
# Residuals 229 4.8148  0.0210    

eta_squared(Chl.lm, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype  |     0.18 | [0.10, 1.00]
# Site      |     0.23 | [0.16, 1.00]
# Origin    | 3.19e-03 | [0.00, 1.00]
# TimeP     |     0.27 | [0.18, 1.00]

emmeans(Chl.lm, pairwise~TimeP)
# contrast estimate     SE  df t.ratio p.value
# M1 - M12  -0.1676 0.0299 229  -5.603  <.0001 *
# M1 - M4    0.0131 0.0298 229   0.442  0.9921
# M1 - M8    0.0508 0.0298 229   1.706  0.4321
# M1 - W2    0.2409 0.0298 229   8.094  <.0001 *
# M12 - M4   0.1808 0.0298 229   6.075  <.0001 *
# M12 - M8   0.2184 0.0298 229   7.340  <.0001 *
# M12 - W2   0.4085 0.0298 229  13.727  <.0001 *
# M4 - M8    0.0376 0.0296 229   1.271  0.7090
# M4 - W2    0.2277 0.0296 229   7.693  <.0001 *
# M8 - W2    0.1901 0.0296 229   6.422  <.0001 *

##Strong Seasonal effect, Split analysis by Timepoints 


#-------Univariate Physiology Week 2------------

Coral.W2<-subset(Coral, TimeP=="W2")

####Protein Host Fraction####

##Check normality
hist(Coral.W2$TP_ug.cm2_C)
qqnorm(Coral.W2$TP_ug.cm2_C)
shapiro.test(Coral.W2$TP_ug.cm2_C)
#Normal

##Model
Prot.C.lm.W2<-lm(TP_ug.cm2_C~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.W2)

##Check residuals
plot(fitted(Prot.C.lm.W2), resid(Prot.C.lm.W2))
abline(0,0)

qqnorm(resid(Prot.C.lm.W2))
qqline(resid(Prot.C.lm.W2))

plot(density(resid(Prot.C.lm.W2)))

##Model results
summary(Prot.C.lm.W2)
# Multiple R-squared:  0.5194,	Adjusted R-squared:  0.4055 
# F-statistic: 4.563 on 9 and 38 DF,  p-value: 0.0004066

anova(Prot.C.lm.W2)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2  58705 29352.5 13.5871 3.536e-05 ***
# Site             1  18812 18812.4  8.7081  0.005401 ** 
# Origin           1   1835  1834.6  0.8492  0.362589    
# Site:Origin      1   2946  2945.5  1.3635  0.250210    
# Genotype:Origin  2   3902  1951.2  0.9032  0.413796    
# Genotype:Site    2   2513  1256.3  0.5815  0.563930    
# Residuals       38  82092  2160.3

eta_squared(Prot.C.lm.W2, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        | 0.34 | [0.13, 1.00]
# Site            | 0.11 | [0.00, 1.00]
# Origin          | 0.01 | [0.00, 1.00]
# Site:Origin     | 0.02 | [0.00, 1.00]
# Genotype:Origin | 0.02 | [0.00, 1.00]
# Genotype:Site   | 0.01 | [0.00, 1.00]

emmeans(Prot.C.lm.W2, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12    -84.6 16.4 38  -5.146  <.0001 *
# AC10 - AC8     -54.1 16.4 38  -3.294  0.0059 *
# AC12 - AC8      30.4 16.4 38   1.852  0.1669


####Protein Symbiont Fraction####

##Check normality
hist(Coral.W2$TP_ug.cm2_S)
qqnorm(Coral.W2$TP_ug.cm2_S)
shapiro.test(Coral.W2$TP_ug.cm2_S)
#Normal

##Model
Prot.S.lm.W2<-lm(TP_ug.cm2_S~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.W2)

##Check residuals
plot(fitted(Prot.S.lm.W2), resid(Prot.S.lm.W2))
abline(0,0)

qqnorm(resid(Prot.S.lm.W2))
qqline(resid(Prot.S.lm.W2))

plot(density(resid(Prot.S.lm.W2)))

##Model results
summary(Prot.S.lm.W2)
# Multiple R-squared:  0.4899,	Adjusted R-squared:  0.3691 
# F-statistic: 4.055 on 9 and 38 DF,  p-value: 0.001049

anova(Prot.S.lm.W2)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2  59802 29901.1  8.0966 0.001177 **
# Site             1     15    15.1  0.0041 0.949395   
# Origin           1  21176 21176.0  5.7340 0.021677 * 
# Site:Origin      1  22463 22462.8  6.0824 0.018279 * 
# Genotype:Origin  2  22121 11060.4  2.9949 0.061973 . 
# Genotype:Site    2   9203  4601.5  1.2460 0.299141   
# Residuals       38 140337  3693.1   

eta_squared(Prot.S.lm.W2, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.22 | [0.04, 1.00]
# Site            | 5.48e-05 | [0.00, 1.00]
# Origin          |     0.08 | [0.00, 1.00]
# Site:Origin     |     0.08 | [0.00, 1.00]
# Genotype:Origin |     0.08 | [0.00, 1.00]
# Genotype:Site   |     0.03 | [0.00, 1.00]

emmeans(Prot.S.lm.W2, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12    -63.8 21.5 38  -2.968  0.0140 *
# AC10 - AC8      18.7 21.5 38   0.870  0.6624
# AC12 - AC8      82.4 21.5 38   3.837  0.0013 *


####Biomass Host Fraction####

##Check normality
hist(Coral.W2$AFDW_mg.cm2_C)
qqnorm(Coral.W2$AFDW_mg.cm2_C)
shapiro.test(Coral.W2$AFDW_mg.cm2_C)
#Normal

##Model
Bio.C.lm.W2<-lm(AFDW_mg.cm2_C~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.W2)

##Check residuals
plot(fitted(Bio.C.lm.W2), resid(Bio.C.lm.W2))
abline(0,0)

qqnorm(resid(Bio.C.lm.W2))
qqline(resid(Bio.C.lm.W2))

plot(density(resid(Bio.C.lm.W2)))

##Model results
summary(Bio.C.lm.W2)
# Multiple R-squared:  0.1975,	Adjusted R-squared:  0.007449 
# F-statistic: 1.039 on 9 and 38 DF,  p-value: 0.4281

anova(Bio.C.lm.W2)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2 0.33781 0.168904  2.2954 0.11452  
# Site             1 0.00119 0.001186  0.0161 0.89965  
# Origin           1 0.31308 0.313084  4.2548 0.04602 *
# Site:Origin      1 0.00683 0.006831  0.0928 0.76227  
# Genotype:Origin  2 0.02783 0.013913  0.1891 0.82850  
# Genotype:Site    2 0.00148 0.000739  0.0100 0.99002  
# Residuals       38 2.79619 0.073584    

eta_squared(Bio.C.lm.W2, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.10 | [0.00, 1.00]
# Site            | 3.40e-04 | [0.00, 1.00]
# Origin          |     0.09 | [0.00, 1.00]
# Site:Origin     | 1.96e-03 | [0.00, 1.00]
# Genotype:Origin | 7.99e-03 | [0.00, 1.00]
# Genotype:Site   | 4.24e-04 | [0.00, 1.00]

emmeans(Bio.C.lm.W2, pairwise~Genotype)
