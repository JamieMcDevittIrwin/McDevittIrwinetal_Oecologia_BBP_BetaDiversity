# Species Richness Linear Models
# CORAL


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
library(sjPlot)
library(influence.ME)
library(lmerTest) # this package extens lmer() to give you p-values
library(broom)

richness.metadata.sc <- read.csv("data/calculated/richness.metadata.scaled.csv")
ls()
#######################################################################


#######################################################################
# Model Selection
######################################################################
# Random Variable Selection 
# when you are comparing models that differ in their random effects use REML

M1 <- lm(richness_coral ~ Net_primary_productivity + Exposure + 
           Scaridae + MG...Bahamas.big.7 + Depth,  
         data = richness.metadata.sc)
M2 <- lmer(richness_coral ~ Net_primary_productivity + Exposure + 
             Scaridae + MG...Bahamas.big.7 + Depth + 
             (1|Reef.name), REML=TRUE, data = richness.metadata.sc)

plot(M1) 

AICc(M1) # 140.4074
AICc(M2) # 141.0339
# Better without the random effect

options(na.action = "na.fail") 
dd <- dredge(M1)
head(dd)
# Model selection table 
#    (Int)    Dpt    Exp Net_prm_prd df  logLik  AICc delta weight
# 11 21.67        1.1430     -0.9372  4 -61.300 132.4  0.00  0.244
# 1  21.67                            2 -64.219 132.9  0.52  0.188
# 3  21.67        0.7894              3 -62.975 133.0  0.58  0.183
# 4  21.67 0.8267 1.0640              4 -61.639 133.1  0.68  0.174
# 12 21.67 0.6403 1.2960     -0.7781  5 -60.462 133.8  1.36  0.124

# best model: exposure + productivity
# next best model: nothing! 
# but its 0.52 AICc apart 

M4 <- lm(richness_coral ~ Net_primary_productivity + Exposure, data = richness.metadata.sc)
summary(M4)
# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               21.6667     0.4783  45.304   <2e-16 ***
# Net_primary_productivity  -0.9372     0.5263  -1.781   0.0876 .  
# Exposure                   1.1433     0.5263   2.172   0.0399 *  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 2.485 on 24 degrees of freedom
# Multiple R-squared:  0.1945,	Adjusted R-squared:  0.1274 
# F-statistic: 2.897 on 2 and 24 DF,  p-value: 0.07463
#######################################################################


#######################################################################
# Model Validation
#######################################################################
# should be done with REML

# 1) Assess Homogeneity (standardized residuals vs fitted values)
plot(M4) 

# 2) QQPlot to Assess Normality 
qqnorm(residuals(M4))
hist(residuals(M4))

# 3) Independence Assumption: plot residuals against each explan var in model (don't want to see a trend, which means they aren't independent)
plot(y = residuals(M4), x = richness.metadata.sc$Net_primary_productivity, xlab = "Net_primary_productivity",ylab = "Residuals") 
abline(0,0)

# 4) Plot residuals against each explan var not used in model 
op <- par(mfrow = c(2, 2))
plot(y = residuals(M4), x = richness.metadata.sc$Depth, xlab = "Depth",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M4), x = richness.metadata.sc$MG...Bahamas.big.7, xlab = "MG...Bahamas.big.7",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M4), x = richness.metadata.sc$Scaridae, xlab = "Scaridae",ylab = "Residuals") 
abline(0,0)
plot(y = residuals(M4), x = richness.metadata.sc$Exposure, xlab = "Exposure",ylab = "Residuals") 
abline(0,0)
par(op)

#######################################################################



#######################################################################
# Do the taxa make the AICc better? 
######################################################################
head(richness.metadata.sc)

# best model
M.test <- lm(richness_coral ~ Net_primary_productivity + Exposure,  data = richness.metadata.sc)
summary(M.test)
AICc(M.test) # 132.4173
coral.div <- tidy(M.test)
save(coral.div , file="results/species-richness/coral-div-summary.Rdata")



predictMODEL = function(DATA, INDEX1, INDEX2, INDEX3){
	# best model with each richness 
	div1 <- lm(richness_coral ~ Net_primary_productivity + Exposure + DATA[, INDEX1], data=DATA)
	div2 <- lm(richness_coral ~ Net_primary_productivity + Exposure + DATA[, INDEX2], data=DATA)
	div3 <- lm(richness_coral ~ Net_primary_productivity + Exposure + DATA[, INDEX3], data=DATA)

	result <- list(AICc(div1), AICc(div2), AICc(div3)) # get AICc of each model 
	names(result) <- c(INDEX1, INDEX2, INDEX3)   # rename
	return(result)
}



results = predictMODEL(richness.metadata.sc, "richness_algae", "richness_sponge", "richness_gorgonian")
results
# not any better 
#######################################################################



