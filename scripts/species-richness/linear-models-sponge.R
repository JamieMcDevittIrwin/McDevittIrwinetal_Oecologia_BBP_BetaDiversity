# Species Richness Linear Models
# SPONGE


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

source("analyses/functions/tidy_lmer.R")

richness.metadata.sc <- read.csv("data/calculated/richness.metadata.scaled.csv")
ls()
#######################################################################



#######################################################################
# UPDATED: May 2020, using poisson for sponge now because of potential heteroskedasticity
#######################################################################

# Trying with GLM and GLMER
M1 <- glm(richness_sponge ~ Net_primary_productivity + Exposure + 
            Scaridae + MG...Bahamas.big.7 + Depth, family=poisson, data = richness.metadata.sc)
M2 <- glmer(richness_sponge ~ Net_primary_productivity + Exposure + 
              Scaridae + MG...Bahamas.big.7 + Depth + (1|Reef.name), family=poisson,
            data = richness.metadata.sc) # singular fit 

AICc(M1) # better model
AICc(M2)

anova(M2, M1) # random effect is not significant
dispersiontest(M1, alternative= "less") # not overdispersed or underdispersed


options(na.action='na.fail')
dd <- dredge(M1)
head(dd)
# depth, MG, prod (1.88 AICc better)


# Best Model
M3 <- glm(richness_sponge ~ Depth + MG...Bahamas.big.7 + Net_primary_productivity,family=poisson, data=richness.metadata.sc)
summary(M3)
AICc(M3) # 175.0435 
dispersiontest(M3) 
library(rsq)
rsq::rsq(M3) # 0.3505798




# Model Validation
plot(M3) 
qqnorm(residuals(M3))
hist(residuals(M3))



# Add other taxonomic groups richness
predictMODEL.glm = function(DATA, INDEX1, INDEX2, INDEX3){
  # best model with each richness 
  div1 <- glm(richness_sponge ~ Depth + MG...Bahamas.big.7 + Net_primary_productivity + DATA[, INDEX1],family=poisson, data=DATA)
  div2 <- glm(richness_sponge ~ Depth + MG...Bahamas.big.7 + Net_primary_productivity + DATA[, INDEX2],family=poisson, data=DATA)
  div3 <- glm(richness_sponge ~ Depth + MG...Bahamas.big.7 + Net_primary_productivity + DATA[, INDEX3],family=poisson, data=DATA)
  
  result <- list(AICc(div1), AICc(div2), AICc(div3)) # get AICc of each model 
  names(result) <- c(INDEX1, INDEX2, INDEX3)   # rename
  return(result)
}


results = predictMODEL.glm(richness.metadata.sc, "richness_algae", "richness_gorgonian", "richness_coral")
results
# improved slightly by gorgonian (1.3 AICc) but not enough to suggest a specific trend
#######################################################################

