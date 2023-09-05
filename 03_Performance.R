#Title: "Coral Performance following Reciprocal Transplant"
#Author: "Serena Hackerott"
#Date: "08/29/2023"

#-------Set Up------------

####Load Packages####
library(ggplot2)
library(effectsize)
library(emmeans)
library(Rmisc)

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


####Biomass Symbiont Fraction####

##Check normality
hist(Coral.W2$AFDW_mg.cm2_S)
qqnorm(Coral.W2$AFDW_mg.cm2_S)
shapiro.test(Coral.W2$AFDW_mg.cm2_S)
#Normal

##Model
Bio.S.lm.W2<-lm(AFDW_mg.cm2_S~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.W2)

##Check residuals
plot(fitted(Bio.S.lm.W2), resid(Bio.S.lm.W2))
abline(0,0)

qqnorm(resid(Bio.S.lm.W2))
qqline(resid(Bio.S.lm.W2))

plot(density(resid(Bio.S.lm.W2)))

##Model results
summary(Bio.S.lm.W2)
# Multiple R-squared:  0.4524,	Adjusted R-squared:  0.3227 
# F-statistic: 3.488 on 9 and 38 DF,  p-value: 0.003153

anova(Bio.S.lm.W2)
# Df  Sum Sq  Mean Sq F value   Pr(>F)   
# Genotype         2 0.13455 0.067274  3.3447 0.045921 * 
# Site             1 0.03109 0.031093  1.5459 0.221364   
# Origin           1 0.12449 0.124491  6.1893 0.017355 * 
# Site:Origin      1 0.15171 0.151707  7.5425 0.009158 **
# Genotype:Origin  2 0.01152 0.005758  0.2863 0.752664   
# Genotype:Site    2 0.17810 0.089051  4.4274 0.018688 * 
# Residuals       38 0.76432 0.020114    

eta_squared(Bio.S.lm.W2, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.10 | [0.00, 1.00]
# Site            |     0.02 | [0.00, 1.00]
# Origin          |     0.09 | [0.00, 1.00]
# Site:Origin     |     0.11 | [0.00, 1.00]
# Genotype:Origin | 8.25e-03 | [0.00, 1.00]
# Genotype:Site   |     0.13 | [0.00, 1.00]


####Chlorophyll####

##Check normality
hist(Coral.W2$Chl_ug.cm2)
qqnorm(Coral.W2$Chl_ug.cm2)
shapiro.test(Coral.W2$Chl_ug.cm2)
#Normal

##Model
Chl.lm.W2<-lm(Chl_ug.cm2~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.W2)

##Check residuals
plot(fitted(Chl.lm.W2), resid(Chl.lm.W2))
abline(0,0)

qqnorm(resid(Chl.lm.W2))
qqline(resid(Chl.lm.W2))

plot(density(resid(Chl.lm.W2)))

##Model results
summary(Chl.lm.W2)
# Multiple R-squared:    0.8,	Adjusted R-squared:  0.7527 
# F-statistic: 16.89 on 9 and 38 DF,  p-value: 9.491e-11

anova(Chl.lm.W2)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 2.85265 1.42632 40.8792 3.373e-10 ***
# Site             1 1.10349 1.10349 31.6267 1.868e-06 ***
# Origin           1 0.32631 0.32631  9.3523  0.004067 ** 
# Site:Origin      1 0.69453 0.69453 19.9056 7.022e-05 ***
# Genotype:Origin  2 0.03502 0.01751  0.5019  0.609329    
# Genotype:Site    2 0.29208 0.14604  4.1856  0.022760 *  
# Residuals       38 1.32586 0.03489       

eta_squared(Chl.lm.W2, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.43 | [0.22, 1.00]
# Site            |     0.17 | [0.03, 1.00]
# Origin          |     0.05 | [0.00, 1.00]
# Site:Origin     |     0.10 | [0.00, 1.00]
# Genotype:Origin | 5.28e-03 | [0.00, 1.00]
# Genotype:Site   |     0.04 | [0.00, 1.00]



#-------Univariate Physiology Month 1------------

Coral.M1<-subset(Coral, TimeP=="M1")

####Protein Host Fraction####

##Check normality
hist(Coral.M1$TP_ug.cm2_C)
qqnorm(Coral.M1$TP_ug.cm2_C)
shapiro.test(Coral.M1$TP_ug.cm2_C)
#Not Normal

hist(log(Coral.M1$TP_ug.cm2_C+1))
qqnorm(log(Coral.M1$TP_ug.cm2_C+1))
shapiro.test(log(Coral.M1$TP_ug.cm2_C+1))
#Not Normal but less skewed

##Model with log(+1) transformation
Prot.C.lm.M1<-lm(log(TP_ug.cm2_C+1)~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M1)

##Check residuals
plot(fitted(Prot.C.lm.M1), resid(Prot.C.lm.M1))
abline(0,0)

qqnorm(resid(Prot.C.lm.M1))
qqline(resid(Prot.C.lm.M1))

plot(density(resid(Prot.C.lm.M1)))

##Model results
summary(Prot.C.lm.M1)
# Multiple R-squared:  0.5224,	Adjusted R-squared:  0.4062 
# F-statistic: 4.496 on 9 and 37 DF,  p-value: 0.0004917

anova(Prot.C.lm.M1)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 0.19941 0.09970  4.8945    0.0130 *  
# Site             1 0.56195 0.56195 27.5862 6.455e-06 ***
# Origin           1 0.00261 0.00261  0.1281    0.7224    
# Site:Origin      1 0.01141 0.01141  0.5600    0.4590    
# Genotype:Origin  2 0.00262 0.00131  0.0642    0.9379    
# Genotype:Site    2 0.04630 0.02315  1.1364    0.3319    
# Residuals       37 0.75372 0.02037   

eta_squared(Prot.C.lm.M1, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.13 | [0.00, 1.00]
# Site            |     0.36 | [0.16, 1.00]
# Origin          | 1.65e-03 | [0.00, 1.00]
# Site:Origin     | 7.23e-03 | [0.00, 1.00]
# Genotype:Origin | 1.66e-03 | [0.00, 1.00]
# Genotype:Site   |     0.03 | [0.00, 1.00]

emmeans(Prot.C.lm.M1, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12  -0.1720 0.0514 37  -3.343  0.0053 *
# AC10 - AC8   -0.0806 0.0505 37  -1.597  0.2596
# AC12 - AC8    0.0914 0.0514 37   1.777  0.1913


####Protein Symbiont Fraction####

##Check normality
hist(Coral.M1$TP_ug.cm2_S)
qqnorm(Coral.M1$TP_ug.cm2_S)
shapiro.test(Coral.M1$TP_ug.cm2_S)
#Normal

##Model
Prot.S.lm.M1<-lm(TP_ug.cm2_S~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M1)

##Check residuals
plot(fitted(Prot.S.lm.M1), resid(Prot.S.lm.M1))
abline(0,0)

qqnorm(resid(Prot.S.lm.M1))
qqline(resid(Prot.S.lm.M1))

plot(density(resid(Prot.S.lm.M1)))

##Model results
summary(Prot.S.lm.M1)
# Multiple R-squared:  0.3805,	Adjusted R-squared:  0.2298 
# F-statistic: 2.525 on 9 and 37 DF,  p-value: 0.02285

anova(Prot.S.lm.M1)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2 105208   52604  7.4981 0.001846 **
# Site             1      4       4  0.0006 0.980422   
# Origin           1     54      54  0.0077 0.930592   
# Site:Origin      1  18803   18803  2.6802 0.110084   
# Genotype:Origin  2  15700    7850  1.1189 0.337424   
# Genotype:Site    2  19670    9835  1.4018 0.258910   
# Residuals       37 259577    7016   

eta_squared(Prot.S.lm.M1, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.25 | [0.06, 1.00]
# Site            | 1.02e-05 | [0.00, 1.00]
# Origin          | 1.29e-04 | [0.00, 1.00]
# Site:Origin     |     0.04 | [0.00, 1.00]
# Genotype:Origin |     0.04 | [0.00, 1.00]
# Genotype:Site   |     0.05 | [0.00, 1.00]

emmeans(Prot.S.lm.M1, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12   -98.21 30.2 37  -3.253  0.0067
# AC10 - AC8      9.96 29.6 37   0.336  0.9397
# AC12 - AC8    108.17 30.2 37   3.583  0.0027


####Biomass Host Fraction####

##Check normality
hist(Coral.M1$AFDW_mg.cm2_C)
qqnorm(Coral.M1$AFDW_mg.cm2_C)
shapiro.test(Coral.M1$AFDW_mg.cm2_C)
#Not Normal

hist(log(Coral.M1$AFDW_mg.cm2_C+1))
qqnorm(log(Coral.M1$AFDW_mg.cm2_C+1))
shapiro.test(log(Coral.M1$AFDW_mg.cm2_C+1))
#Still not Normal, but less skewed

##Model with log(+1) transformation
Bio.C.lm.M1<-lm(log(AFDW_mg.cm2_C+1)~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M1)

##Check residuals
plot(fitted(Bio.C.lm.M1), resid(Bio.C.lm.M1))
abline(0,0)

qqnorm(resid(Bio.C.lm.M1))
qqline(resid(Bio.C.lm.M1))

plot(density(resid(Bio.C.lm.M1)))

##Model results
summary(Bio.C.lm.M1)
# Multiple R-squared:  0.4666,	Adjusted R-squared:  0.3369 
# F-statistic: 3.597 on 9 and 37 DF,  p-value: 0.002665

anova(Bio.C.lm.M1)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2 0.18546 0.092732  8.7683 0.0007639 ***
# Site             1 0.03003 0.030028  2.8393 0.1004002    
# Origin           1 0.00469 0.004693  0.4438 0.5094315    
# Site:Origin      1 0.00239 0.002386  0.2256 0.6375884    
# Genotype:Origin  2 0.00143 0.000715  0.0676 0.9347762    
# Genotype:Site    2 0.11834 0.059171  5.5949 0.0075345 ** 
# Residuals       37 0.39131 0.010576    

eta_squared(Bio.C.lm.M1, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.25 | [0.06, 1.00]
# Site            |     0.04 | [0.00, 1.00]
# Origin          | 6.40e-03 | [0.00, 1.00]
# Site:Origin     | 3.25e-03 | [0.00, 1.00]
# Genotype:Origin | 1.95e-03 | [0.00, 1.00]
# Genotype:Site   |     0.16 | [0.01, 1.00]


####Biomass Symbiont Fraction####

##Check normality
hist(Coral.M1$AFDW_mg.cm2_S)
qqnorm(Coral.M1$AFDW_mg.cm2_S)
shapiro.test(Coral.M1$AFDW_mg.cm2_S)
#Normal

##Model
Bio.S.lm.M1<-lm(AFDW_mg.cm2_S~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M1)

##Check residuals
plot(fitted(Bio.S.lm.M1), resid(Bio.S.lm.M1))
abline(0,0)

qqnorm(resid(Bio.S.lm.M1))
qqline(resid(Bio.S.lm.M1))

plot(density(resid(Bio.S.lm.M1)))

##Model results
summary(Bio.S.lm.M1)
# Multiple R-squared:  0.1129,	Adjusted R-squared:  -0.1029 
# F-statistic: 0.5233 on 9 and 37 DF,  p-value: 0.8481

anova(Bio.S.lm.M1)
# Df  Sum Sq  Mean Sq F value   Pr(>F)   
# Genotype         2 0.01175 0.005875  0.1552 0.8568
# Site             1 0.04177 0.041771  1.1033 0.3004
# Origin           1 0.03824 0.038236  1.0099 0.3215
# Site:Origin      1 0.00490 0.004900  0.1294 0.7211
# Genotype:Origin  2 0.02813 0.014065  0.3715 0.6923
# Genotype:Site    2 0.05353 0.026765  0.7069 0.4997
# Residuals       37 1.40089 0.037862     

eta_squared(Bio.S.lm.M1, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        | 7.44e-03 | [0.00, 1.00]
# Site            |     0.03 | [0.00, 1.00]
# Origin          |     0.02 | [0.00, 1.00]
# Site:Origin     | 3.10e-03 | [0.00, 1.00]
# Genotype:Origin |     0.02 | [0.00, 1.00]
# Genotype:Site   |     0.03 | [0.00, 1.00]


####Chlorophyll####

##Check normality
hist(Coral.M1$Chl_ug.cm2)
qqnorm(Coral.M1$Chl_ug.cm2)
shapiro.test(Coral.M1$Chl_ug.cm2)
#Not Normal

hist(log(Coral.M1$Chl_ug.cm2+1))
qqnorm(log(Coral.M1$Chl_ug.cm2+1))
shapiro.test(log(Coral.M1$Chl_ug.cm2+1))
#Normal

##Model with log(+1) transformation
Chl.lm.M1<-lm(log(Chl_ug.cm2+1)~Genotype+Site+Origin+
                Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M1)

##Check residuals
plot(fitted(Chl.lm.M1), resid(Chl.lm.M1))
abline(0,0)

qqnorm(resid(Chl.lm.M1))
qqline(resid(Chl.lm.M1))

plot(density(resid(Chl.lm.M1)))

##Model results
summary(Chl.lm.M1)
# Multiple R-squared:  0.6543,	Adjusted R-squared:  0.5702 
# F-statistic: 7.782 on 9 and 37 DF,  p-value: 2.54e-06

anova(Chl.lm.M1)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 0.52180 0.260902 17.2708 5.039e-06 ***
# Site             1 0.20060 0.200605 13.2794 0.0008186 ***
# Origin           1 0.01137 0.011369  0.7526 0.3912584    
# Site:Origin      1 0.00752 0.007522  0.4979 0.4848335    
# Genotype:Origin  2 0.02452 0.012258  0.8114 0.4519661    
# Genotype:Site    2 0.29218 0.146089  9.6706 0.0004183 ***
# Residuals       37 0.55894 0.015107      

eta_squared(Chl.lm.M1, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.32 | [0.11, 1.00]
# Site            |     0.12 | [0.01, 1.00]
# Origin          | 7.03e-03 | [0.00, 1.00]
# Site:Origin     | 4.65e-03 | [0.00, 1.00]
# Genotype:Origin |     0.02 | [0.00, 1.00]
# Genotype:Site   |     0.18 | [0.01, 1.00]



#-------Univariate Physiology Month 4------------

Coral.M4<-subset(Coral, TimeP=="M4")

####Protein Host Fraction####

##Check normality
hist(Coral.M4$TP_ug.cm2_C)
qqnorm(Coral.M4$TP_ug.cm2_C)
shapiro.test(Coral.M4$TP_ug.cm2_C)
#Normal

##Model 
Prot.C.lm.M4<-lm(TP_ug.cm2_C~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M4)

##Check residuals
plot(fitted(Prot.C.lm.M4), resid(Prot.C.lm.M4))
abline(0,0)

qqnorm(resid(Prot.C.lm.M4))
qqline(resid(Prot.C.lm.M4))

plot(density(resid(Prot.C.lm.M4)))

##Model results
summary(Prot.C.lm.M4)
# Multiple R-squared:  0.4026,	Adjusted R-squared:  0.2611 
# F-statistic: 2.845 on 9 and 38 DF,  p-value: 0.01155

anova(Prot.C.lm.M4)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2  15990  7994.9  4.4158 0.018864 * 
# Site             1  22904 22903.8 12.6505 0.001026 **
# Origin           1      3     3.5  0.0019 0.965167   
# Site:Origin      1    614   614.3  0.3393 0.563691   
# Genotype:Origin  2   4334  2166.9  1.1969 0.313273   
# Genotype:Site    2   2512  1256.0  0.6937 0.505931   
# Residuals       38  68799  1810.5     

eta_squared(Prot.C.lm.M4, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.14 | [0.00, 1.00]
# Site            |     0.20 | [0.04, 1.00]
# Origin          | 3.04e-05 | [0.00, 1.00]
# Site:Origin     | 5.33e-03 | [0.00, 1.00]
# Genotype:Origin |     0.04 | [0.00, 1.00]
# Genotype:Site   |     0.02 | [0.00, 1.00]

emmeans(Prot.C.lm.M4, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12    -44.7 15 38  -2.970  0.0139 *
# AC10 - AC8     -21.2 15 38  -1.408  0.3470
# AC12 - AC8      23.5 15 38   1.562  0.2741

####Protein Symbiont Fraction####

##Check normality
hist(Coral.M4$TP_ug.cm2_S)
qqnorm(Coral.M4$TP_ug.cm2_S)
shapiro.test(Coral.M4$TP_ug.cm2_S)
#Normal

##Model
Prot.S.lm.M4<-lm(TP_ug.cm2_S~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M4)

##Check residuals
plot(fitted(Prot.S.lm.M4), resid(Prot.S.lm.M4))
abline(0,0)

qqnorm(resid(Prot.S.lm.M4))
qqline(resid(Prot.S.lm.M4))

plot(density(resid(Prot.S.lm.M4)))

##Model results
summary(Prot.S.lm.M4)
# Multiple R-squared:  0.4681,	Adjusted R-squared:  0.3421 
# F-statistic: 3.716 on 9 and 38 DF,  p-value: 0.002016

anova(Prot.S.lm.M4)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2  31872   15936  5.2709 0.0095432 ** 
# Site             1  44406   44406 14.6875 0.0004627 ***
# Origin           1   1696    1696  0.5610 0.4584847    
# Site:Origin      1   3711    3711  1.2276 0.2748410    
# Genotype:Origin  2   4584    2292  0.7580 0.4755344    
# Genotype:Site    2  14847    7424  2.4554 0.0993400 .  
# Residuals       38 114889    3023 

eta_squared(Prot.S.lm.M4, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.15 | [0.00, 1.00]
# Site            |     0.21 | [0.05, 1.00]
# Origin          | 7.85e-03 | [0.00, 1.00]
# Site:Origin     |     0.02 | [0.00, 1.00]
# Genotype:Origin |     0.02 | [0.00, 1.00]
# Genotype:Site   |     0.07 | [0.00, 1.00]


emmeans(Prot.S.lm.M4, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12      -10 19.4 38  -0.514  0.8649
# AC10 - AC8        49 19.4 38   2.519  0.0416 *
# AC12 - AC8        59 19.4 38   3.033  0.0118 *


####Biomass Host Fraction####

##Check normality
hist(Coral.M4$AFDW_mg.cm2_C)
qqnorm(Coral.M4$AFDW_mg.cm2_C)
shapiro.test(Coral.M4$AFDW_mg.cm2_C)
#Normal

##Model
Bio.C.lm.M4<-lm(AFDW_mg.cm2_C~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M4)

##Check residuals
plot(fitted(Bio.C.lm.M4), resid(Bio.C.lm.M4))
abline(0,0)

qqnorm(resid(Bio.C.lm.M4))
qqline(resid(Bio.C.lm.M4))

plot(density(resid(Bio.C.lm.M4)))

##Model results
summary(Bio.C.lm.M4)
# Multiple R-squared:  0.254,	Adjusted R-squared:  0.07737 
# F-statistic: 1.438 on 9 and 38 DF,  p-value: 0.2069

anova(Bio.C.lm.M4)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2 0.13516 0.067579  2.2078 0.12385  
# Site             1 0.00442 0.004425  0.1446 0.70591  
# Origin           1 0.00973 0.009725  0.3177 0.57630  
# Site:Origin      1 0.20377 0.203775  6.6573 0.01386 *
# Genotype:Origin  2 0.04268 0.021340  0.6972 0.50424  
# Genotype:Site    2 0.00037 0.000183  0.0060 0.99405  
# Residuals       38 1.16315 0.030609     

eta_squared(Bio.C.lm.M4, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.09 | [0.00, 1.00]
# Site            | 2.84e-03 | [0.00, 1.00]
# Origin          | 6.24e-03 | [0.00, 1.00]
# Site:Origin     |     0.13 | [0.01, 1.00]
# Genotype:Origin |     0.03 | [0.00, 1.00]
# Genotype:Site   | 2.34e-04 | [0.00, 1.00]


####Biomass Symbiont Fraction####

##Check normality
hist(Coral.M4$AFDW_mg.cm2_S)
qqnorm(Coral.M4$AFDW_mg.cm2_S)
shapiro.test(Coral.M4$AFDW_mg.cm2_S)
#Normal

##Model
Bio.S.lm.M4<-lm(AFDW_mg.cm2_S~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M4)

##Check residuals
plot(fitted(Bio.S.lm.M4), resid(Bio.S.lm.M4))
abline(0,0)

qqnorm(resid(Bio.S.lm.M4))
qqline(resid(Bio.S.lm.M4))

plot(density(resid(Bio.S.lm.M4)))

##Model results
summary(Bio.S.lm.M4)
# Multiple R-squared:  0.5427,	Adjusted R-squared:  0.4344 
# F-statistic: 5.012 on 9 and 38 DF,  p-value: 0.0001812

anova(Bio.S.lm.M4)
# Df  Sum Sq  Mean Sq F value   Pr(>F)   
# Genotype         2 0.04176 0.02088  1.5079   0.23432    
# Site             1 0.40470 0.40470 29.2251 3.713e-06 ***
# Origin           1 0.05124 0.05124  3.7004   0.06191 .  
# Site:Origin      1 0.03643 0.03643  2.6310   0.11307    
# Genotype:Origin  2 0.00756 0.00378  0.2730   0.76259    
# Genotype:Site    2 0.08289 0.04145  2.9931   0.06207 .  
# Residuals       38 0.52621 0.01385       

eta_squared(Bio.S.lm.M4, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.04 | [0.00, 1.00]
# Site            |     0.35 | [0.16, 1.00]
# Origin          |     0.04 | [0.00, 1.00]
# Site:Origin     |     0.03 | [0.00, 1.00]
# Genotype:Origin | 6.57e-03 | [0.00, 1.00]
# Genotype:Site   |     0.07 | [0.00, 1.00]


####Chlorophyll####

##Check normality
hist(Coral.M4$Chl_ug.cm2)
qqnorm(Coral.M4$Chl_ug.cm2)
shapiro.test(Coral.M4$Chl_ug.cm2)
#Not Normal

hist(log(Coral.M4$Chl_ug.cm2+1))
qqnorm(log(Coral.M4$Chl_ug.cm2+1))
shapiro.test(log(Coral.M4$Chl_ug.cm2+1))
#Normal

##Model with log(+1) transformation
Chl.lm.M4<-lm(log(Chl_ug.cm2+1)~Genotype+Site+Origin+
                Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M4)

##Check residuals
plot(fitted(Chl.lm.M4), resid(Chl.lm.M4))
abline(0,0)

qqnorm(resid(Chl.lm.M4))
qqline(resid(Chl.lm.M4))

plot(density(resid(Chl.lm.M4)))

##Model results
summary(Chl.lm.M4)
# Multiple R-squared:  0.7726,	Adjusted R-squared:  0.7188 
# F-statistic: 14.35 on 9 and 38 DF,  p-value: 9.718e-10

anova(Chl.lm.M4)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 0.11987 0.05993  4.0975  0.024467 *  
# Site             1 1.34967 1.34967 92.2732 1.031e-11 ***
# Origin           1 0.00008 0.00008  0.0057  0.940319    
# Site:Origin      1 0.16733 0.16733 11.4398  0.001678 ** 
# Genotype:Origin  2 0.03021 0.01510  1.0327  0.365829    
# Genotype:Site    2 0.22139 0.11069  7.5679  0.001712 ** 
# Residuals       38 0.55582 0.01463       

eta_squared(Chl.lm.M4, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.05 | [0.00, 1.00]
# Site            |     0.55 | [0.37, 1.00]
# Origin          | 3.40e-05 | [0.00, 1.00]
# Site:Origin     |     0.07 | [0.00, 1.00]
# Genotype:Origin |     0.01 | [0.00, 1.00]
# Genotype:Site   |     0.09 | [0.00, 1.00]


#-------Univariate Physiology Month 8------------

Coral.M8<-subset(Coral, TimeP=="M8")

####Protein Host Fraction####

##Check normality
hist(Coral.M8$TP_ug.cm2_C)
qqnorm(Coral.M8$TP_ug.cm2_C)
shapiro.test(Coral.M8$TP_ug.cm2_C)
#Normal

##Model 
Prot.C.lm.M8<-lm(TP_ug.cm2_C~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M8)

##Check residuals
plot(fitted(Prot.C.lm.M8), resid(Prot.C.lm.M8))
abline(0,0)

qqnorm(resid(Prot.C.lm.M8))
qqline(resid(Prot.C.lm.M8))

plot(density(resid(Prot.C.lm.M8)))

##Model results
summary(Prot.C.lm.M8)
# Multiple R-squared:  0.3472,	Adjusted R-squared:  0.1926 
# F-statistic: 2.246 on 9 and 38 DF,  p-value: 0.03994

anova(Prot.C.lm.M8)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2  65154   32577  7.3682 0.001976 **
# Site             1   2649    2649  0.5991 0.443695   
# Origin           1   7495    7495  1.6951 0.200770   
# Site:Origin      1    134     134  0.0304 0.862489   
# Genotype:Origin  2   8103    4052  0.9164 0.408625   
# Genotype:Site    2   5819    2910  0.6581 0.523632   
# Residuals       38 168011    4421     

eta_squared(Prot.C.lm.M8, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.25 | [0.06, 1.00]
# Site            |     0.01 | [0.00, 1.00]
# Origin          |     0.03 | [0.00, 1.00]
# Site:Origin     | 5.22e-04 | [0.00, 1.00]
# Genotype:Origin |     0.03 | [0.00, 1.00]
# Genotype:Site   |     0.02 | [0.00, 1.00]

emmeans(Prot.C.lm.M8, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12    -87.2 23.5 38  -3.710  0.0019 *
# AC10 - AC8     -63.7 23.5 38  -2.708  0.0267 *
# AC12 - AC8      23.6 23.5 38   1.003  0.5798

####Protein Symbiont Fraction####

##Check normality
hist(Coral.M8$TP_ug.cm2_S)
qqnorm(Coral.M8$TP_ug.cm2_S)
shapiro.test(Coral.M8$TP_ug.cm2_S)
#Normal

##Model
Prot.S.lm.M8<-lm(TP_ug.cm2_S~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M8)

##Check residuals
plot(fitted(Prot.S.lm.M8), resid(Prot.S.lm.M8))
abline(0,0)

qqnorm(resid(Prot.S.lm.M8))
qqline(resid(Prot.S.lm.M8))

plot(density(resid(Prot.S.lm.M8)))

##Model results
summary(Prot.S.lm.M8)
# Multiple R-squared:  0.4294,	Adjusted R-squared:  0.2943 
# F-statistic: 3.177 on 9 and 38 DF,  p-value: 0.00587

anova(Prot.S.lm.M8)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2  26586   13293  3.5389 0.0389571 *  
#   Site             1  66510   66510 17.7070 0.0001516 ***
#   Origin           1   1415    1415  0.3766 0.5430809    
# Site:Origin      1   4903    4903  1.3052 0.2604119    
# Genotype:Origin  2    906     453  0.1206 0.8867697    
# Genotype:Site    2   7094    3547  0.9443 0.3978714    
# Residuals       38 142734    3756 

eta_squared(Prot.S.lm.M8, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.11 | [0.00, 1.00]
# Site            |     0.27 | [0.09, 1.00]
# Origin          | 5.66e-03 | [0.00, 1.00]
# Site:Origin     |     0.02 | [0.00, 1.00]
# Genotype:Origin | 3.62e-03 | [0.00, 1.00]
# Genotype:Site   |     0.03 | [0.00, 1.00]

emmeans(Prot.S.lm.M8, pairwise~Genotype)
# contrast    estimate   SE df t.ratio p.value
# AC10 - AC12    -57.1 21.7 38  -2.637  0.0316 *
# AC10 - AC8     -35.2 21.7 38  -1.626  0.2471
# AC12 - AC8      21.9 21.7 38   1.010  0.5751


####Biomass Host Fraction####

##Check normality
hist(Coral.M8$AFDW_mg.cm2_C)
qqnorm(Coral.M8$AFDW_mg.cm2_C)
shapiro.test(Coral.M8$AFDW_mg.cm2_C)
#Normal

##Model
Bio.C.lm.M8<-lm(AFDW_mg.cm2_C~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M8)

##Check residuals
plot(fitted(Bio.C.lm.M8), resid(Bio.C.lm.M8))
abline(0,0)

qqnorm(resid(Bio.C.lm.M8))
qqline(resid(Bio.C.lm.M8))

plot(density(resid(Bio.C.lm.M8)))

##Model results
summary(Bio.C.lm.M8)
# Multiple R-squared:  0.3022,	Adjusted R-squared:  0.137 
# F-statistic: 1.829 on 9 and 38 DF,  p-value: 0.0946

anova(Bio.C.lm.M8)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2 0.32133 0.160666  3.3215 0.04683 *
# Site             1 0.29141 0.291413  6.0246 0.01880 *
# Origin           1 0.04683 0.046826  0.9681 0.33139  
# Site:Origin      1 0.04384 0.043840  0.9063 0.34711  
# Genotype:Origin  2 0.07314 0.036569  0.7560 0.47647  
# Genotype:Site    2 0.01960 0.009798  0.2026 0.81751  
# Residuals       38 1.83809 0.048371  

eta_squared(Bio.C.lm.M8, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.12 | [0.00, 1.00]
# Site            |     0.11 | [0.00, 1.00]
# Origin          |     0.02 | [0.00, 1.00]
# Site:Origin     |     0.02 | [0.00, 1.00]
# Genotype:Origin |     0.03 | [0.00, 1.00]
# Genotype:Site   | 7.44e-03 | [0.00, 1.00]


####Biomass Symbiont Fraction####

##Check normality
hist(Coral.M8$AFDW_mg.cm2_S)
qqnorm(Coral.M8$AFDW_mg.cm2_S)
shapiro.test(Coral.M8$AFDW_mg.cm2_S)
#Not Normal

hist(log(Coral.M8$AFDW_mg.cm2_S+1))
qqnorm(log(Coral.M8$AFDW_mg.cm2_S+1))
shapiro.test(log(Coral.M8$AFDW_mg.cm2_S+1))
#Normal

##Model with log(+1) transformation
Bio.S.lm.M8<-lm(log(AFDW_mg.cm2_S+1)~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M8)

##Check residuals
plot(fitted(Bio.S.lm.M8), resid(Bio.S.lm.M8))
abline(0,0)

qqnorm(resid(Bio.S.lm.M8))
qqline(resid(Bio.S.lm.M8))

plot(density(resid(Bio.S.lm.M8)))

##Model results
summary(Bio.S.lm.M8)
# Multiple R-squared:  0.3411,	Adjusted R-squared:  0.185 
# F-statistic: 2.185 on 9 and 38 DF,  p-value: 0.04527

anova(Bio.S.lm.M8)
# Df  Sum Sq  Mean Sq F value   Pr(>F)   
# Genotype         2 0.07764 0.038818  3.3289 0.046540 * 
# Site             1 0.11315 0.113150  9.7034 0.003491 **
# Origin           1 0.00003 0.000030  0.0026 0.959652   
# Site:Origin      1 0.00493 0.004927  0.4225 0.519595   
# Genotype:Origin  2 0.02266 0.011332  0.9718 0.387620   
# Genotype:Site    2 0.01095 0.005473  0.4693 0.629003   
# Residuals       38 0.44311 0.011661            

eta_squared(Bio.S.lm.M8, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.12 | [0.00, 1.00]
# Site            |     0.17 | [0.03, 1.00]
# Origin          | 4.50e-05 | [0.00, 1.00]
# Site:Origin     | 7.33e-03 | [0.00, 1.00]
# Genotype:Origin |     0.03 | [0.00, 1.00]
# Genotype:Site   |     0.02 | [0.00, 1.00]


####Chlorophyll####

##Check normality
hist(Coral.M8$Chl_ug.cm2)
qqnorm(Coral.M8$Chl_ug.cm2)
shapiro.test(Coral.M8$Chl_ug.cm2)
#Not Normal

hist(log(Coral.M8$Chl_ug.cm2+1))
qqnorm(log(Coral.M8$Chl_ug.cm2+1))
shapiro.test(log(Coral.M8$Chl_ug.cm2+1))
#Normal

##Model with log(+1) transformation
Chl.lm.M8<-lm(log(Chl_ug.cm2+1)~Genotype+Site+Origin+
                Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M8)

##Check residuals
plot(fitted(Chl.lm.M8), resid(Chl.lm.M8))
abline(0,0)

qqnorm(resid(Chl.lm.M8))
qqline(resid(Chl.lm.M8))

plot(density(resid(Chl.lm.M8)))

##Model results
summary(Chl.lm.M8)
# Multiple R-squared:  0.7163,	Adjusted R-squared:  0.6491 
# F-statistic: 10.66 on 9 and 38 DF,  p-value: 5.082e-08

anova(Chl.lm.M8)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 0.66085 0.33042 19.7183 1.336e-06 ***
# Site             1 0.74263 0.74263 44.3171 7.184e-08 ***
# Origin           1 0.01593 0.01593  0.9507   0.33571    
# Site:Origin      1 0.05601 0.05601  3.3427   0.07536 .  
# Genotype:Origin  2 0.04070 0.02035  1.2143   0.30819    
# Genotype:Site    2 0.09168 0.04584  2.7355   0.07764 .  
# Residuals       38 0.63677 0.01676       

eta_squared(Chl.lm.M8, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.05 | [0.00, 1.00]
# Site            |     0.55 | [0.37, 1.00]
# Origin          | 3.40e-05 | [0.00, 1.00]
# Site:Origin     |     0.07 | [0.00, 1.00]
# Genotype:Origin |     0.01 | [0.00, 1.00]
# Genotype:Site   |     0.09 | [0.00, 1.00]



#-------Univariate Physiology Month 12------------

Coral.M12<-subset(Coral, TimeP=="M12")

####Protein Host Fraction####

##Check normality
hist(Coral.M12$TP_ug.cm2_C)
qqnorm(Coral.M12$TP_ug.cm2_C)
shapiro.test(Coral.M12$TP_ug.cm2_C)
#Normal

##Model 
Prot.C.lm.M12<-lm(TP_ug.cm2_C~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M12)

##Check residuals
plot(fitted(Prot.C.lm.M12), resid(Prot.C.lm.M12))
abline(0,0)

qqnorm(resid(Prot.C.lm.M12))
qqline(resid(Prot.C.lm.M12))

plot(density(resid(Prot.C.lm.M12)))

##Model results
summary(Prot.C.lm.M12)
# Multiple R-squared:  0.1394,	Adjusted R-squared:  -0.07571 
# F-statistic: 0.6481 on 9 and 36 DF,  p-value: 0.7486

anova(Prot.C.lm.M12)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2  11623  5811.6  0.9338 0.4024
# Site             1   3830  3830.1  0.6154 0.4379
# Origin           1   1955  1955.5  0.3142 0.5786
# Site:Origin      1  13000 13000.0  2.0888 0.1570
# Genotype:Origin  2   4829  2414.3  0.3879 0.6813
# Genotype:Site    2   1064   532.0  0.0855 0.9183
# Residuals       36 224050  6223.6    


####Protein Symbiont Fraction####

##Check normality
hist(Coral.M12$TP_ug.cm2_S)
qqnorm(Coral.M12$TP_ug.cm2_S)
shapiro.test(Coral.M12$TP_ug.cm2_S)
#Normal

##Model
Prot.S.lm.M12<-lm(TP_ug.cm2_S~Genotype+Site+Origin+
                   Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M12)

##Check residuals
plot(fitted(Prot.S.lm.M12), resid(Prot.S.lm.M12))
abline(0,0)

qqnorm(resid(Prot.S.lm.M12))
qqline(resid(Prot.S.lm.M12))

plot(density(resid(Prot.S.lm.M12)))

##Model results
summary(Prot.S.lm.M12)
# Multiple R-squared:  0.5267,	Adjusted R-squared:  0.4115 
# F-statistic: 4.574 on 9 and 37 DF,  p-value: 0.0004266

anova(Prot.S.lm.M12)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2  17574    8787  2.8948   0.06792 .  
# Site             1  97874   97874 32.2431 1.713e-06 ***
#   Origin           1    989     989  0.3258   0.57158    
# Site:Origin      1    116     116  0.0381   0.84635    
# Genotype:Origin  2   1462     731  0.2408   0.78725    
# Genotype:Site    2   6957    3479  1.1460   0.32893    
# Residuals       37 112314    3036 

eta_squared(Prot.S.lm.M12, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.07 | [0.00, 1.00]
# Site            |     0.41 | [0.21, 1.00]
# Origin          | 4.17e-03 | [0.00, 1.00]
# Site:Origin     | 4.87e-04 | [0.00, 1.00]
# Genotype:Origin | 6.16e-03 | [0.00, 1.00]
# Genotype:Site   |     0.03 | [0.00, 1.00]


####Biomass Host Fraction####

##Check normality
hist(Coral.M12$AFDW_mg.cm2_C)
qqnorm(Coral.M12$AFDW_mg.cm2_C)
shapiro.test(Coral.M12$AFDW_mg.cm2_C)
#Normal

##Model
Bio.C.lm.M12<-lm(AFDW_mg.cm2_C~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M12)

##Check residuals
plot(fitted(Bio.C.lm.M12), resid(Bio.C.lm.M12))
abline(0,0)

qqnorm(resid(Bio.C.lm.M12))
qqline(resid(Bio.C.lm.M12))

plot(density(resid(Bio.C.lm.M12)))

##Model results
summary(Bio.C.lm.M12)
# Multiple R-squared:  0.3103,	Adjusted R-squared:  0.1425 
# F-statistic:  1.85 on 9 and 37 DF,  p-value: 0.09159

anova(Bio.C.lm.M12)
# Df Sum Sq Mean Sq F value   Pr(>F)   
# Genotype         2 0.2208 0.11038  1.1882 0.31614  
# Site             1 0.3321 0.33211  3.5748 0.06651 .
# Origin           1 0.0138 0.01378  0.1483 0.70233  
# Site:Origin      1 0.3622 0.36219  3.8986 0.05583 .
# Genotype:Origin  2 0.4338 0.21690  2.3347 0.11095  
# Genotype:Site    2 0.1839 0.09195  0.9898 0.38129  
# Residuals       37 3.4375 0.09290    


####Biomass Symbiont Fraction####

##Check normality
hist(Coral.M12$AFDW_mg.cm2_S)
qqnorm(Coral.M12$AFDW_mg.cm2_S)
shapiro.test(Coral.M12$AFDW_mg.cm2_S)
#Normal

##Model
Bio.S.lm.M12<-lm(AFDW_mg.cm2_S~Genotype+Site+Origin+
                  Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M12)

##Check residuals
plot(fitted(Bio.S.lm.M12), resid(Bio.S.lm.M12))
abline(0,0)

qqnorm(resid(Bio.S.lm.M12))
qqline(resid(Bio.S.lm.M12))

plot(density(resid(Bio.S.lm.M12)))

##Model results
summary(Bio.S.lm.M12)
# Multiple R-squared:  0.1743,	Adjusted R-squared:  -0.02657 
# F-statistic: 0.8677 on 9 and 37 DF,  p-value: 0.5616

anova(Bio.S.lm.M12)
# Df  Sum Sq  Mean Sq F value   Pr(>F)   
# Genotype         2 0.01081 0.005403  0.1356 0.87366  
# Site             1 0.25175 0.251748  6.3157 0.01645 *
# Origin           1 0.00533 0.005330  0.1337 0.71670  
# Site:Origin      1 0.00206 0.002056  0.0516 0.82160  
# Genotype:Origin  2 0.01007 0.005034  0.1263 0.88175  
# Genotype:Site    2 0.03129 0.015643  0.3924 0.67818  
# Residuals       37 1.47484 0.039861              

eta_squared(Bio.S.lm.M12, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        | 6.05e-03 | [0.00, 1.00]
# Site            |     0.14 | [0.01, 1.00]
# Origin          | 2.98e-03 | [0.00, 1.00]
# Site:Origin     | 1.15e-03 | [0.00, 1.00]
# Genotype:Origin | 5.64e-03 | [0.00, 1.00]
# Genotype:Site   |     0.02 | [0.00, 1.00]


####Chlorophyll####

##Check normality
hist(Coral.M12$Chl_ug.cm2)
qqnorm(Coral.M12$Chl_ug.cm2)
shapiro.test(Coral.M12$Chl_ug.cm2)
#Not Normal

hist(log(Coral.M12$Chl_ug.cm2+1))
qqnorm(log(Coral.M12$Chl_ug.cm2+1))
shapiro.test(log(Coral.M12$Chl_ug.cm2+1))
#Normal

##Model with log(+1) transformation
Chl.lm.M12<-lm(log(Chl_ug.cm2+1)~Genotype+Site+Origin+
                Site:Origin + Genotype:Origin+ Genotype:Site, data=Coral.M12)

##Check residuals
plot(fitted(Chl.lm.M12), resid(Chl.lm.M12))
abline(0,0)

qqnorm(resid(Chl.lm.M12))
qqline(resid(Chl.lm.M12))

plot(density(resid(Chl.lm.M12)))

##Model results
summary(Chl.lm.M12)
# Multiple R-squared:  0.8144,	Adjusted R-squared:  0.7692 
# F-statistic: 18.04 on 9 and 37 DF,  p-value: 5.215e-11

anova(Chl.lm.M12)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Genotype         2 0.90347 0.45174 27.7442 4.356e-08 ***
# Site             1 1.59963 1.59963 98.2440 5.838e-12 ***
# Origin           1 0.00080 0.00080  0.0494   0.82531    
# Site:Origin      1 0.01398 0.01398  0.8588   0.36009    
# Genotype:Origin  2 0.02598 0.01299  0.7978   0.45792    
# Genotype:Site    2 0.09940 0.04970  3.0525   0.05929 .  
# Residuals       37 0.60244 0.01628        

eta_squared(Chl.lm.M12, partial=FALSE)
# Parameter |     Eta2 |       95% CI
# Genotype        |     0.28 | [0.08, 1.00]
# Site            |     0.49 | [0.30, 1.00]
# Origin          | 2.48e-04 | [0.00, 1.00]
# Site:Origin     | 4.31e-03 | [0.00, 1.00]
# Genotype:Origin | 8.00e-03 | [0.00, 1.00]
# Genotype:Site   |     0.03 | [0.00, 1.00]

####Summary####

##Week 2
#Protein of Symbiont, Biomass of Symbiont, and Chlorophyll: Significant Origin and Site: Origin Effects

##Month 1
##Only Genotype, Site, and GenoXSite effects

##Month 4 
##SitexOrigin effect for Biomass of Host and Chlorophyll, marginal (p<0.1) for Biomass of Symbiont

##Month 8 
##Marginal (p<0.1) SitexOrigin effect for Chlorophyll

##Month 12
##Only Genotype, Site, and GenoXSite effects


####Figures- Chlorophyll####

Geno.Site.colors.o<-c("#466BE3FF",   "#28BBECFF",  "#1AE4B6FF", 
                                 "#A2FC3CFF","#9D1001FF", "#F26014FF")
                  

# ("#466BE3FF", "#466BE3FF", "#28BBECFF", "#28BBECFF", 
# "#1AE4B6FF", "#1AE4B6FF", "#A2FC3CFF", "#A2FC3CFF",
# "#9D1001FF", "#9D1001FF", "#F26014FF", "#F26014FF")
                                   
####Week 2

##Summary Statistics 
Chl.SE.W2<-summarySE(Coral.W2, measurevar="Chl_ug.cm2", groupvars=c("Genotype", "Site", "Orig", "Origin"), na.rm=TRUE)

str(Chl.SE.W2)
Chl.SE.W2$Site<-factor(Chl.SE.W2$Site, levels=c("SS", "KL"), ordered=TRUE)
Chl.SE.W2$Genotype<-factor(Chl.SE.W2$Genotype, levels=c("AC8", "AC10", "AC12"), ordered=TRUE)
Chl.SE.W2$Orig<-factor(Chl.SE.W2$Orig, levels=c("N", "T"), ordered=TRUE)
Chl.SE.W2$Origin<-factor(Chl.SE.W2$Origin, levels=c("Native", "Transplant"), ordered=TRUE)

Chl.SE.W2$Set<-paste(Chl.SE.W2$Site, Chl.SE.W2$Genotype, Chl.SE.W2$Orig, sep=".")
Chl.SE.W2$Geno.Orig<-paste(Chl.SE.W2$Genotype, Chl.SE.W2$Orig, sep=".")
Chl.SE.W2$Geno.Site<-paste(Chl.SE.W2$Genotype, Chl.SE.W2$Site,  sep=".")
Chl.SE.W2$Geno.Site<-factor(Chl.SE.W2$Geno.Site, levels=c("AC8.SS", "AC8.KL", "AC10.SS",  "AC10.KL", "AC12.SS",  "AC12.KL"), ordered=TRUE)

##Plot Mean +/- Standard Deviation by Genotype and Site
Chl.SE.W2.plot<-ggplot(Chl.SE.W2, aes(x=Geno.Site, y=Chl_ug.cm2, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_ug.cm2-sd, ymax=Chl_ug.cm2+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 4.1)+
  labs(x="Genotype and Site", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")));Chl.SE.W2.plot


####Month 1

##Summary Statistics 
Chl.SE.M1<-summarySE(Coral.M1, measurevar="Chl_ug.cm2", groupvars=c("Genotype", "Site", "Orig", "Origin"), na.rm=TRUE)

str(Chl.SE.M1)
Chl.SE.M1$Site<-factor(Chl.SE.M1$Site, levels=c("SS", "KL"), ordered=TRUE)
Chl.SE.M1$Genotype<-factor(Chl.SE.M1$Genotype, levels=c("AC8", "AC10", "AC12"), ordered=TRUE)
Chl.SE.M1$Orig<-factor(Chl.SE.M1$Orig, levels=c("N", "T"), ordered=TRUE)
Chl.SE.M1$Origin<-factor(Chl.SE.M1$Origin, levels=c("Native", "Transplant"), ordered=TRUE)

Chl.SE.M1$Set<-paste(Chl.SE.M1$Site, Chl.SE.M1$Genotype, Chl.SE.M1$Orig, sep=".")
Chl.SE.M1$Geno.Orig<-paste(Chl.SE.M1$Genotype, Chl.SE.M1$Orig, sep=".")
Chl.SE.M1$Geno.Site<-paste(Chl.SE.M1$Genotype, Chl.SE.M1$Site,  sep=".")
Chl.SE.M1$Geno.Site<-factor(Chl.SE.M1$Geno.Site, levels=c("AC8.SS", "AC8.KL", "AC10.SS",  "AC10.KL", "AC12.SS",  "AC12.KL"), ordered=TRUE)

##Plot Mean +/- Standard Deviation by Genotype and Site
Chl.SE.M1.plot<-ggplot(Chl.SE.M1, aes(x=Geno.Site, y=Chl_ug.cm2, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_ug.cm2-sd, ymax=Chl_ug.cm2+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 4.1)+
  labs(x="Genotype and Site", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")));Chl.SE.M1.plot


####Month 4

##Summary Statistics 
Chl.SE.M4<-summarySE(Coral.M4, measurevar="Chl_ug.cm2", groupvars=c("Genotype", "Site", "Orig", "Origin"), na.rm=TRUE)

str(Chl.SE.M4)
Chl.SE.M4$Site<-factor(Chl.SE.M4$Site, levels=c("SS", "KL"), ordered=TRUE)
Chl.SE.M4$Genotype<-factor(Chl.SE.M4$Genotype, levels=c("AC8", "AC10", "AC12"), ordered=TRUE)
Chl.SE.M4$Orig<-factor(Chl.SE.M4$Orig, levels=c("N", "T"), ordered=TRUE)
Chl.SE.M4$Origin<-factor(Chl.SE.M4$Origin, levels=c("Native", "Transplant"), ordered=TRUE)

Chl.SE.M4$Set<-paste(Chl.SE.M4$Site, Chl.SE.M4$Genotype, Chl.SE.M4$Orig, sep=".")
Chl.SE.M4$Geno.Orig<-paste(Chl.SE.M4$Genotype, Chl.SE.M4$Orig, sep=".")
Chl.SE.M4$Geno.Site<-paste(Chl.SE.M4$Genotype, Chl.SE.M4$Site,  sep=".")
Chl.SE.M4$Geno.Site<-factor(Chl.SE.M4$Geno.Site, levels=c("AC8.SS", "AC8.KL", "AC10.SS",  "AC10.KL", "AC12.SS",  "AC12.KL"), ordered=TRUE)

##Plot Mean +/- Standard Deviation by Genotype and Site
Chl.SE.M4.plot<-ggplot(Chl.SE.M4, aes(x=Geno.Site, y=Chl_ug.cm2, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_ug.cm2-sd, ymax=Chl_ug.cm2+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 4.1)+
  labs(x="Genotype and Site", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")));Chl.SE.M4.plot



####Month 8

##Summary Statistics 
Chl.SE.M8<-summarySE(Coral.M8, measurevar="Chl_ug.cm2", groupvars=c("Genotype", "Site", "Orig", "Origin"), na.rm=TRUE)

str(Chl.SE.M8)
Chl.SE.M8$Site<-factor(Chl.SE.M8$Site, levels=c("SS", "KL"), ordered=TRUE)
Chl.SE.M8$Genotype<-factor(Chl.SE.M8$Genotype, levels=c("AC8", "AC10", "AC12"), ordered=TRUE)
Chl.SE.M8$Orig<-factor(Chl.SE.M8$Orig, levels=c("N", "T"), ordered=TRUE)
Chl.SE.M8$Origin<-factor(Chl.SE.M8$Origin, levels=c("Native", "Transplant"), ordered=TRUE)

Chl.SE.M8$Set<-paste(Chl.SE.M8$Site, Chl.SE.M8$Genotype, Chl.SE.M8$Orig, sep=".")
Chl.SE.M8$Geno.Orig<-paste(Chl.SE.M8$Genotype, Chl.SE.M8$Orig, sep=".")
Chl.SE.M8$Geno.Site<-paste(Chl.SE.M8$Genotype, Chl.SE.M8$Site,  sep=".")
Chl.SE.M8$Geno.Site<-factor(Chl.SE.M8$Geno.Site, levels=c("AC8.SS", "AC8.KL", "AC10.SS",  "AC10.KL", "AC12.SS",  "AC12.KL"), ordered=TRUE)

##Plot Mean +/- Standard Deviation by Genotype and Site
Chl.SE.M8.plot<-ggplot(Chl.SE.M8, aes(x=Geno.Site, y=Chl_ug.cm2, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_ug.cm2-sd, ymax=Chl_ug.cm2+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 4.1)+
  labs(x="Genotype and Site", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")));Chl.SE.M8.plot


####Month 12

##Summary Statistics 
Chl.SE.M12<-summarySE(Coral.M12, measurevar="Chl_ug.cm2", groupvars=c("Genotype", "Site", "Orig", "Origin"), na.rm=TRUE)

str(Chl.SE.M12)
Chl.SE.M12$Site<-factor(Chl.SE.M12$Site, levels=c("SS", "KL"), ordered=TRUE)
Chl.SE.M12$Genotype<-factor(Chl.SE.M12$Genotype, levels=c("AC8", "AC10", "AC12"), ordered=TRUE)
Chl.SE.M12$Orig<-factor(Chl.SE.M12$Orig, levels=c("N", "T"), ordered=TRUE)
Chl.SE.M12$Origin<-factor(Chl.SE.M12$Origin, levels=c("Native", "Transplant"), ordered=TRUE)

Chl.SE.M12$Set<-paste(Chl.SE.M12$Site, Chl.SE.M12$Genotype, Chl.SE.M12$Orig, sep=".")
Chl.SE.M12$Geno.Orig<-paste(Chl.SE.M12$Genotype, Chl.SE.M12$Orig, sep=".")
Chl.SE.M12$Geno.Site<-paste(Chl.SE.M12$Genotype, Chl.SE.M12$Site,  sep=".")
Chl.SE.M12$Geno.Site<-factor(Chl.SE.M12$Geno.Site, levels=c("AC8.SS", "AC8.KL", "AC10.SS",  "AC10.KL", "AC12.SS",  "AC12.KL"), ordered=TRUE)

##Plot Mean +/- Standard Deviation by Genotype and Site
Chl.SE.M12.plot<-ggplot(Chl.SE.M12, aes(x=Geno.Site, y=Chl_ug.cm2, colour=Geno.Site, group=Origin)) + 
  geom_errorbar(aes(ymin=Chl_ug.cm2-sd, ymax=Chl_ug.cm2+sd), width=cap.sz, position=position_dodge(0.3), size=bar.sz) +
  geom_point(aes(shape=Origin), position=position_dodge(0.3), size=point.sz, stroke=bar.sz+1)+
  theme_classic()+
  scale_colour_manual(values=Geno.Site.colors.o)+
  scale_shape_manual(values=c(19, 1))+
  theme(axis.title.x = element_text(size = axis.title.sz), axis.title.y = element_text(size = axis.title.sz), 
        axis.text.x=element_text(size=axis.txt.sz, colour="black"), axis.text.y=element_text(size=axis.txt.sz, colour="black"), 
        legend.text=element_text(size=leg.txt.sz), legend.title=element_text(size=leg.title.sz), 
        legend.box.background = element_rect(color = "black"))+
  ylim(0, 4.1)+
  labs(x="Genotype and Site", y=expression(paste('Chlorophyll-a (\u03BCg cm'^-2*")")));Chl.SE.M12.plot

