# MRM models with Sorensen beta diversity #
# ***************************************************
# ndividual env and linear distance to find significant drivers 
# then adding other functional groups to see if they are significant 
# ***************************************************
setwd("/Users/jamiemcdevitt-irwin/Documents/Documents - Jamie's MacBook Pro/Git_Repos/Stanford_BBP")
getwd()

# clear my environment
rm(list=ls()) 
library(stringr)
library(dplyr)
library(tidyr)
library(geosphere)
library(vegan)
library(tibble)
library(ecodist)
library(ggplot2)
library(egg)

load("data/metadata/site_linear_distance.Rdata")

load("data/calculated/scaridae_disimilarity_matrix.Rdata")
load("data/calculated/productivity_disimilarity_matrix.Rdata")
load("data/calculated/exposure_disimilarity_matrix.Rdata")
load("data/calculated/marketgravity_disimilarity_matrix.Rdata")
load("data/calculated/depth_disimilarity_matrix.Rdata")

load("data/calculated/montastraea_stonycoral_sorensen_matrix.Rdata")
load("data/calculated/montastraea_macroalgae_sorensen_matrix.Rdata")
load("data/calculated/montastraea_sponge_sorensen_matrix.Rdata")
load("data/calculated/montastraea_gorgonian_sorensen_matrix.Rdata")
ls()

###############################################################


###############################################################
# MRM Function 

# Coral
coral_mrm1 <- MRM(stony.coral.sor ~ depth_distance + exposure_dist + mg7_dist + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
coral_mrm1 # sig predictors= depth

coral_mrm2 <- MRM(stony.coral.sor ~ depth_distance, nperm=9999)
coral_mrm2 # all sig

# Pull out the sig ones for rmarkdown results
coral_coefs_sor <-data.frame(coral_mrm2[1])
coral_R_sor <- data.frame(coral_mrm2[2])
coral_F_sor <- data.frame(coral_mrm2[3])

# Algae
algae_mrm1 <- MRM(macroalgae.sor ~ depth_distance + exposure_dist 
	+ mg7_dist + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
algae_mrm1  # significant predictors: depth, exposure, productivity

algae_mrm2 <- MRM(macroalgae.sor ~ depth_distance + exposure_dist + productivity_dist, nperm=9999)
algae_mrm2  # significant predictors

# Pull out the sig ones for rmarkdown results
algae_coefs_sor <-data.frame(algae_mrm2[1])
algae_R_sor <- data.frame(algae_mrm2[2])
algae_F_sor <- data.frame(algae_mrm2[3])

# Sponge
sponge_mrm1 <- MRM(sponge.sor ~ depth_distance + exposure_dist + mg7_dist + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
sponge_mrm1 # significant predictors: none

sponge_mrm2 <- MRM(sponge.sor ~ 1, nperm=9999)

# Gorgonian
gorg_mrm1 <- MRM(gorgonian.sor ~ depth_distance + exposure_dist + mg7_dist + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
gorg_mrm1 # significant predictors: exposure (almost depth)

gorg_mrm2 <- MRM(gorgonian.sor ~ exposure_dist, nperm=9999)
gorg_mrm2

# Pull out the sig ones for rmarkdown results
gorg_coefs_sor <-data.frame(gorg_mrm2[1])
gorg_R_sor <- data.frame(gorg_mrm2[2])
gorg_F_sor <- data.frame(gorg_mrm2[3])


# Save results for rmarkdown read in
ind.mrm.coefs.sor <- list(coral_coefs_sor,algae_coefs_sor,gorg_coefs_sor)
ind.mrm.R.sor <- list(coral_R_sor,algae_R_sor,gorg_R_sor)
ind.mrm.F.sor <- list(coral_F_sor,algae_F_sor,gorg_F_sor)

save(ind.mrm.coefs.sor, file="results/beta-diversity/individual-mrm-environment-results-coefs-sor.Rdata")
save(ind.mrm.R.sor, file="results/beta-diversity/individual-mrm-environment-results-R-sor.Rdata")
save(ind.mrm.F.sor, file="results/beta-diversity/individual-mrm-environment-results-F-sor.Rdata")
###############################################################

###############################################################
# Best model with other functional group added

# Coral
coral_mrm_best <- MRM(stony.coral.sor ~ depth_distance, nperm=9999)
coral_mrm_algae <- MRM(stony.coral.sor ~ depth_distance + macroalgae.sor, nperm=9999)
coral_mrm_sponge <- MRM(stony.coral.sor ~ depth_distance + sponge.sor, nperm=9999)
coral_mrm_gorg <- MRM(stony.coral.sor ~ depth_distance + gorgonian.sor, nperm=9999)

coral_mrm_best # 
coral_mrm_algae # macroalgae is not significatn
coral_mrm_sponge # sponge is not significant
coral_mrm_gorg # gorgonian is not significant


# Algae
algae_mrm_best <- MRM(macroalgae.sor ~ depth_distance + exposure_dist + productivity_dist, nperm=9999)
algae_mrm_coral <- MRM(macroalgae.sor ~ depth_distance + exposure_dist + productivity_dist + stony.coral.sor, nperm=9999)
algae_mrm_sponge <- MRM(macroalgae.sor ~ depth_distance + exposure_dist + productivity_dist + sponge.sor, nperm=9999)
algae_mrm_gorg <- MRM(macroalgae.sor ~ depth_distance + exposure_dist + productivity_dist + gorgonian.sor, nperm=9999)

algae_mrm_best # 0.38
algae_mrm_coral # coral is not significant
algae_mrm_sponge #  sponge is not significant 
algae_mrm_gorg #  gorgonian is not significant


# Sponge
# best model has no predictors
sponge_mrm_coral <- MRM(sponge.sor ~ stony.coral.sor, nperm=9999)
sponge_mrm_algae <- MRM(sponge.sor ~ macroalgae.sor, nperm=9999)
sponge_mrm_gorg <- MRM(sponge.sor ~ gorgonian.sor, nperm=9999)

sponge_mrm_coral # 0.19, coral is not significant
sponge_mrm_algae # 0.25, algae is not significant
sponge_mrm_gorg # gorgonian is not signficant 


# Gorgonian
gorg_mrm_best <- MRM(gorgonian.sor ~  exposure_dist, nperm=9999)
gorg_mrm_coral <- MRM(gorgonian.sor ~  exposure_dist + stony.coral.sor, nperm=9999)
gorg_mrm_algae <- MRM(gorgonian.sor ~ exposure_dist + macroalgae.sor, nperm=9999)
gorg_mrm_sponge <- MRM(gorgonian.sor ~ exposure_dist + sponge.sor, nperm=9999)


gorg_mrm_best
gorg_mrm_coral # coral is not significant
gorg_mrm_algae # algae is not significant
gorg_mrm_sponge # sponge is not significant
###############################################################



###############################################################
# Pull out the coefficients and p-values of the best env model to plot

# Convert to dataframe 
coefDATA = function(DATA){
	ha <- as.data.frame(DATA$coef)
	ha <- rownames_to_column(ha)
	colnames(ha) <- c("env", "coef", "pval")
	return(ha)

}

# Coral
coral.plot = coefDATA(coral_mrm_best)
group <- c("coral", "coral")
coral.plot$group <- group

# Algae
algae.plot = coefDATA(algae_mrm_best)
group <- c("algae", "algae", "algae", "algae")
algae.plot$group <- group

# Sponge has no significant environmental predictors 
# # Sponge
# sponge.plot = coefDATA(sponge_mrm_best)
# group <- c("sponge", "sponge")
# sponge.plot$group <- group

# Gorgonian
gorg.plot = coefDATA(gorg_mrm_best)
group <- c("gorgonian", "gorgonian")
gorg.plot$group <- group

coefs.env <- rbind(coral.plot, algae.plot, gorg.plot)
coefs.env

# add column for label
coral.plot$label <- "e) Coral"
algae.plot$label <- "f) Algal"
gorg.plot$label <- "g) Gorgonian"

coefs.env <- rbind(coral.plot, algae.plot, gorg.plot)
head(coefs.env)

coefs.env <- coefs.env %>% filter(env != "Int")


# change order for the plot
coefs.env$group = factor(coefs.env$group, levels=c('coral','algae','gorgonian'))

###############################################################


###############################################################
# PLOT! 

library(scales)
show_col(hue_pal()(6))


p= ggplot(coefs.env, aes(x=env, y=coef, color=env)) + 
  geom_point(size=8) + 
  facet_grid(group ~ .) + 
  theme(
	panel.grid.major = element_blank(),
	panel.border = element_rect(colour = "black", fill= NA), 
	panel.grid.minor = element_blank(), 
		panel.background= element_blank(), 
	axis.text=element_text(size=11), 
	axis.title=element_text(size=13), 
		axis.text.y=element_blank(),
	strip.text = element_text(size = 11),
	#legend.position= "none",
	plot.title=element_text(face= "bold")) + 
  labs(colour = "Environmental Beta Diversity", 
       x= "", 
       y= "  ",
       subtitle= "Presence-Absence") +
	coord_flip() + geom_hline(yintercept = 0, linetype=2) +
  ylim(-0.005,0.17) +
	scale_colour_manual(values=c("#F8766D", "#B79F00", "#619CFF"),
		breaks=c("depth_distance", "exposure_dist", "productivity_dist"),
		labels=c("Depth", "Exposure", "Productivity"))+  
  ggtitle("Environment")

my_tag <- c("(i) Coral", "(j) Algal", "(k) Gorgonian")
p = tag_facet(p, tag_pool = my_tag,
              open = "", close = "",
              fontface =1,
              hjust = c(-.89, -.88, -0.5),x=Inf, y=-Inf)
p

ggsave(file='figures/montastraea/beta-diversity/mrm-env-coefs-sorensen.pdf', 
       plot=p,
       height=7, width=8)

# save to combine with other figure
mrm_sorenson <- p
mrm_sorenson
save(mrm_sorenson, file='data/figures/mrm_sorenson_env.Rdata')
###############################################################





