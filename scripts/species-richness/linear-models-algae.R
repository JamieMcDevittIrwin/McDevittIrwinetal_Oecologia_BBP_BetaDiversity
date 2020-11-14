# Species Richness Linear Models
# ALGAE 

setwd("/Users/jamiemcdevitt-irwin/Documents/Documents - Jamie's MacBook Pro/Git_Repos/Stanford_BBP")
getwd()

# clear my environment
rm(list=ls()) 
library(stringr)
library(dplyr)
library(tidyr)
library(vegan)
library(lme4)
library(tibble)
library(MuMIn)
library(MASS)
library(influence.ME)
library(broom)
library(AER)

# Purpose: Run linear models on the four functional group species richness 

richness.metadata.sc <- read.csv("data/calculated/richness.metadata.scaled.csv")
ls()
#######################################################################



#######################################################################
# Model Selection
######################################################################
# Random Variable Selection 
# when you are comparing models that differ in their random effects use REML
M1 <- lm(richness_algae ~ Net_primary_productivity + Exposure + 
           Scaridae + MG...Bahamas.big.7 + Depth,  data = richness.metadata.sc)
M2 <- lmer(richness_algae ~ Net_primary_productivity + Exposure + 
             Scaridae + MG...Bahamas.big.7 + Depth + (1|Reef.name), REML=TRUE, 
           data = richness.metadata.sc)

AICc(M1) # 117.4052
AICc(M2) # 122.6937
# Better without the random effect 

# plot residuals vs fitted for full model
plot(M1) 


# Model Selection
options(na.action=na.fail)
dd <- dredge(M1)  
subset(dd, delta < 4)
# Model selection table 
#    (Int)     Dpt     Exp MG...Bhm.big.7 Net_prm_prd      Scr df  logLik  AICc
# 1  17.67                                                      2 -52.282 109.1
# 5  17.67                         0.3371                       3 -51.747 110.5
# 2  17.67 -0.3144                                              3 -51.818 110.7
# 23 17.67         -0.9468         1.3620             -0.77820  5 -48.963 110.8
# 9  17.67                                     0.2057           3 -52.086 111.2
# 7  17.67         -0.5637         0.6907                       4 -50.788 111.4
# 3  17.67         -0.1304                                      3 -52.204 111.5
# 17 17.67                                            -0.07001  3 -52.260 111.6
# 21 17.67                         0.5423             -0.37040  4 -51.282 112.4
# 4  17.67 -0.4022 -0.2641                                      4 -51.518 112.9
# 6  17.67 -0.2111                 0.2498                       4 -51.568 113.0
#    delta weight
# 1   0.00  0.245
# 5   1.47  0.117
# 2   1.61  0.109
# 23  1.72  0.104


# Best model: just the intercept 
M3 <- lm(richness_algae ~ 1, data=richness.metadata.sc)
summary(M3) 

#######################################################################



#######################################################################
# Model Validation
#######################################################################
# should be done with REML

# 1) Assess Homogeneity (standardized residuals vs fitted values)
plot(M3) 

# 2) QQPlot to Assess Normality 
qqnorm(residuals(M3))
hist(residuals(M3))

# 3) Independence Assumption: plot residuals against each explan var in model (don't want to see a trend, which means they aren't independent)

# 3) Plot residuals against each explan var not used in model 
op <- par(mfrow = c(3, 2))
plot(y = residuals(M3), x = richness.metadata.sc$MG...Bahamas.big.7, xlab = "MG...Bahamas.big.7",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M3), x = richness.metadata.sc$Scaridae, xlab = "Scaridae",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M3), x = richness.metadata.sc$Exposure, xlab = "Exposure",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M3), x = richness.metadata.sc$Depth, xlab = "Depth",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M3), x = richness.metadata.sc$Net_primary_productivity, xlab = "Net_primary_productivity",ylab = "Residuals") 
abline(0,0)
par(op)

#######################################################################



#######################################################################
# Add in other functional group richness
######################################################################
head(richness.metadata.sc)

# best model
M3 <- lm(richness_algae ~ 1, data=richness.metadata.sc)
AICc(M3) # 109.0649
algae.div <- tidy(M3)
save(algae.div , file="results/species-richness/algae-div-summary.Rdata")


predictMODEL = function(DATA, INDEX1, INDEX2, INDEX3){
	# best model with each richness 
	div1 <- lm(richness_algae ~ DATA[, INDEX1], data=DATA)
	div2 <- lm(richness_algae ~ DATA[, INDEX2], data=DATA)
	div3 <- lm(richness_algae ~ DATA[, INDEX3], data=DATA)

	result <- list(AICc(div1), AICc(div2), AICc(div3)) # get AICc of each model 
	names(result) <- c(INDEX1, INDEX2, INDEX3)   # rename
	return(result)
}



results = predictMODEL(richness.metadata.sc, "richness_coral", "richness_sponge", "richness_gorgonian")
results

# model is not improved by any of them 
#######################################################################




