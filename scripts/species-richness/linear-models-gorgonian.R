# Species Richness Linear Models
# GORGONIAN

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
library(lmerTest) # this package extens lmer() to give you p-values
library(broom)

# Purpose: Run linear models on the four functional group species richness 

richness.metadata.sc <- read.csv("data/calculated/richness.metadata.scaled.csv")
ls()
#######################################################################



#######################################################################
# Model Selection
######################################################################
# Random Variable Selection 
# when you are comparing models that differ in their random effects use REML

M1 <- lm(richness_gorgonian ~ Net_primary_productivity + Exposure + 
           Scaridae + MG...Bahamas.big.7 + Depth,  
         data = richness.metadata.sc)
M2 <- lmer(richness_gorgonian ~ Net_primary_productivity + Exposure + 
             Scaridae + MG...Bahamas.big.7 + Depth + 
             (1|Reef.name), REML=TRUE, data = richness.metadata.sc)

plot(M1)
plot(M2)

AICc(M1) # 103.2369
AICc(M2) # 111.8757
# Better without the random effect  

# Model Selection
options(na.action="na.fail")
dredge(M1) 

# Model selection table 
#    (Int)     Dpt   Exp MG...Bhm.big.7 Net_prm_prd      Scr df  logLik  AICc
# 31 7.074         3.021        -1.4020      0.6656  0.78070  6 -42.585 101.4
# 32 7.074 -0.3457 3.014        -1.5380      0.6005  0.82090  7 -41.671 103.2
# 15 7.074         2.654        -0.7105      0.5904           5 -45.383 103.6


# best model: Net_primary_productivity + Exposure + Scaridae + MG...Bahamas.big.7
M3 <- lm(richness_gorgonian ~ Net_primary_productivity + Exposure + 
           Scaridae + MG...Bahamas.big.7,  data = richness.metadata.sc)
summary(M3) # Adjusted R-squared:  0.7933 


#######################################################################


#######################################################################
# Model Validation
#######################################################################
# should be done with REML (if you have a random effect)

# 1) Assess Homogeneity (standardized residuals vs fitted values)
plot(M3) 

# 2) QQPlot to Assess Normality 
qqnorm(residuals(M3))
hist(residuals(M3)) 

# 3) Independence Assumption: plot residuals against each explan var in model (don't want to see a trend, which means they aren't independent)
op <- par(mfrow = c(2, 2))
plot(y = residuals(M3), x = richness.metadata.sc$MG...Bahamas.big.7, xlab = "MG...Bahamas.big.7",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M3), x = richness.metadata.sc$Scaridae, xlab = "Scaridae",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M3), x = richness.metadata.sc$Exposure, xlab = "Exposure",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M3), x = richness.metadata.sc$Net_primary_productivity, xlab = "Net_primary_productivity",ylab = "Residuals") 
abline(0,0)
par(op)



# 3) Plot residuals against each explan var not used in model 
plot(y = residuals(M3), x = richness.metadata.sc$Depth, xlab = "Depth",ylab = "Residuals") 
abline(0,0)

#######################################################################



#######################################################################
# Add in other functional group richness
######################################################################
head(richness.metadata.sc)

# best model
M3 <- lm(richness_gorgonian ~ Net_primary_productivity + Exposure + 
           Scaridae + MG...Bahamas.big.7,  data = richness.metadata.sc)
AICc(M3) # 101.3696
gorg.div <- tidy(M3)
save(gorg.div , file="results/species-richness/gorgonian-div-summary.Rdata")



predictMODEL = function(DATA, INDEX1, INDEX2, INDEX3){
	# best model with each richness 
	div1 <- lm(richness_gorgonian ~ Net_primary_productivity + Exposure + Scaridae + MG...Bahamas.big.7 + DATA[, INDEX1], data=DATA)
	div2 <- lm(richness_gorgonian ~ Net_primary_productivity + Exposure + Scaridae + MG...Bahamas.big.7 + DATA[, INDEX2], data=DATA)
	div3 <- lm(richness_gorgonian ~ Net_primary_productivity + Exposure + Scaridae + MG...Bahamas.big.7 + DATA[, INDEX3], data=DATA)

	result <- list(AICc(div1), AICc(div2), AICc(div3)) # get AICc of each model 
	names(result) <- c(INDEX1, INDEX2, INDEX3)   # rename
	return(result)
}



results = predictMODEL(richness.metadata.sc, "richness_algae", "richness_sponge", "richness_coral")
results

# model is not improved by any of them 
# model is slightly improved by algae (by 0.09 AICc) but this isn't large enough
#######################################################################









