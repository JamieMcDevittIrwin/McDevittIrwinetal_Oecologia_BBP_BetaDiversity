# MRM models with Bray Curtis # 
# ***************************************************
# individual env and linear distance to find significant drivers 
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
library(broom)
library(egg)

load("data/metadata/site_linear_distance.Rdata")

load("data/calculated/scaridae_disimilarity_matrix.Rdata")
load("data/calculated/productivity_disimilarity_matrix.Rdata")
load("data/calculated/exposure_disimilarity_matrix.Rdata")
load("data/calculated/marketgravity_disimilarity_matrix.Rdata")
load("data/calculated/depth_disimilarity_matrix.Rdata")

load("data/calculated/montastraea_stonycoral_bray_matrix.Rdata")
load("data/calculated/montastraea_macroalgae_bray_matrix.Rdata")
load("data/calculated/montastraea_sponge_bray_matrix.Rdata")
load("data/calculated/montastraea_gorgonian_bray_matrix.Rdata")
ls()

###############################################################


###############################################################
# MRM Function 

# Coral
coral_mrm1 <- MRM(stony.coral.bray ~ depth_distance + 
                    exposure_dist + mg7_dist + productivity_dist + 
                    scaridae_dist + linear_dist, nperm=9999)
coral_mrm1 # depth, productivity, scaridae, linear distance


coral_mrm2 <- MRM(stony.coral.bray ~ depth_distance + 
                    productivity_dist + 
                    scaridae_dist + 
                    linear_dist, nperm=9999)
coral_mrm2 # all sig

coral_coefs <-data.frame(coral_mrm2[1])
coral_R <- data.frame(coral_mrm2[2])
coral_F <- data.frame(coral_mrm2[3])

# Algae
algae_mrm1 <- MRM(macroalgae.bray ~ depth_distance + exposure_dist + mg7_dist + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
algae_mrm1  # significant predictors: depth (0.05), MG, exposure

algae_mrm2 <- MRM(macroalgae.bray ~ depth_distance + mg7_dist + exposure_dist, nperm=9999)
algae_mrm2  # significant predictors

algae_coefs <-data.frame(algae_mrm2[1])
algae_R <- data.frame(algae_mrm2[2])
algae_F <- data.frame(algae_mrm2[3])

# Sponge
sponge_mrm1 <- MRM(sponge.bray ~ depth_distance + exposure_dist + mg7_dist + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
sponge_mrm1 # significant predictors: mg7 

sponge_mrm2 <- MRM(sponge.bray ~ mg7_dist, nperm=9999)
sponge_mrm2 # significant predictors: mg7 

sponge_coefs <-data.frame(sponge_mrm2[1])
sponge_R <- data.frame(sponge_mrm2[2])
sponge_F <- data.frame(sponge_mrm2[3])

# Gorgonian
gorg_mrm1 <- MRM(gorgonian.bray ~ depth_distance + exposure_dist + mg7_dist + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
gorg_mrm1 # significant predictors: depth, exposure,linear distance , mg7

gorg_mrm2 <- MRM(gorgonian.bray ~ depth_distance + exposure_dist + linear_dist + mg7_dist, nperm=9999)
gorg_mrm2

gorg_coefs <-data.frame(gorg_mrm2[1])
gorg_R <- data.frame(gorg_mrm2[2])
gorg_F <- data.frame(gorg_mrm2[3])

# Save results for rmarkdown read in
ind.mrm.coefs <- list(coral_coefs,algae_coefs,sponge_coefs,gorg_coefs)
ind.mrm.R <- list(coral_R,algae_R,sponge_R,gorg_R)
ind.mrm.F <- list(coral_F,algae_F,sponge_F,gorg_F)

save(ind.mrm.coefs, file="results/beta-diversity/individual-mrm-environment-results-coefs.Rdata")
save(ind.mrm.R, file="results/beta-diversity/individual-mrm-environment-results-R.Rdata")
save(ind.mrm.F, file="results/beta-diversity/individual-mrm-environment-results-F.Rdata")
###############################################################



###############################################################
# Best model with other functional group added
###############################################################

### Coral ###
coral_mrm_best <- MRM(stony.coral.bray ~ depth_distance + productivity_dist + scaridae_dist + linear_dist, nperm=9999)
coral_mrm_algae <- MRM(stony.coral.bray ~ depth_distance + productivity_dist + scaridae_dist + linear_dist + macroalgae.bray, nperm=9999)
coral_mrm_sponge <- MRM(stony.coral.bray ~ depth_distance + productivity_dist + scaridae_dist + linear_dist + sponge.bray, nperm=9999)
coral_mrm_gorg <- MRM(stony.coral.bray ~ depth_distance + productivity_dist + scaridae_dist + linear_dist + gorgonian.bray, nperm=9999)

coral_mrm_best # 0.45
coral_mrm_algae # 0.46, macroalgae is significant
coral_mrm_sponge # 0.46, sponge is significant
coral_mrm_gorg # gorgonian is not significant

# Pull out the sig ones for rmarkdown results
coral_mrm_algae_coefs <-data.frame(coral_mrm_algae[1])[6,]
coral_mrm_algae_R <- data.frame(coral_mrm_algae[2])
coral_mrm_algae_F <- data.frame(coral_mrm_algae[3])

coral_mrm_sponge_coefs <-data.frame(coral_mrm_sponge[1])[6,]
coral_mrm_sponge_R <- data.frame(coral_mrm_sponge[2])
coral_mrm_sponge_F <- data.frame(coral_mrm_sponge[3])



### Algae ###
algae_mrm_best <- MRM(macroalgae.bray ~ depth_distance + mg7_dist+ exposure_dist, nperm=9999)
algae_mrm_coral <- MRM(macroalgae.bray ~ depth_distance + mg7_dist + exposure_dist+ stony.coral.bray, nperm=9999)
algae_mrm_sponge <- MRM(macroalgae.bray ~ depth_distance + mg7_dist + exposure_dist+ sponge.bray, nperm=9999)
algae_mrm_gorg <- MRM(macroalgae.bray ~ depth_distance + mg7_dist + exposure_dist+ gorgonian.bray, nperm=9999)

algae_mrm_best # 0.21
algae_mrm_coral # 0.26, coral is significant
algae_mrm_sponge # 0.30, sponge is significant 
algae_mrm_gorg # 0.31, gorgonian is significant

# Pull out the sig ones for rmarkdown results
algae_mrm_coral_coefs <-data.frame(algae_mrm_coral[1])[5,]
algae_mrm_coral_R <- data.frame(algae_mrm_coral[2])
algae_mrm_coral_F <- data.frame(algae_mrm_coral[3])

algae_mrm_sponge_coefs <-data.frame(algae_mrm_sponge[1])[5,]
algae_mrm_sponge_R <- data.frame(algae_mrm_sponge[2])
algae_mrm_sponge_F <- data.frame(algae_mrm_sponge[3])

algae_mrm_gorg_coefs <-data.frame(algae_mrm_gorg[1])[5,]
algae_mrm_gorg_R <- data.frame(algae_mrm_gorg[2])
algae_mrm_gorg_F <- data.frame(algae_mrm_gorg[3])



### Sponge ###
sponge_mrm_best <- MRM(sponge.bray ~ mg7_dist, nperm=9999)
sponge_mrm_coral <- MRM(sponge.bray ~ mg7_dist + stony.coral.bray, nperm=9999)
sponge_mrm_algae <- MRM(sponge.bray ~ mg7_dist + macroalgae.bray, nperm=9999)
sponge_mrm_gorg <- MRM(sponge.bray ~ mg7_dist + gorgonian.bray, nperm=9999)

sponge_mrm_best # 0.15
sponge_mrm_coral # 0.19, coral is significant
sponge_mrm_algae # 0.25, algae is significant
sponge_mrm_gorg # gorgonian is not signficant 

# Pull out the sig ones for rmarkdown results
sponge_mrm_coral_coefs <-data.frame(sponge_mrm_coral[1])[3,]
sponge_mrm_coral_R <- data.frame(sponge_mrm_coral[2])
sponge_mrm_coral_F <- data.frame(sponge_mrm_coral[3])

sponge_mrm_algae_coefs <-data.frame(sponge_mrm_algae[1])[3,]
sponge_mrm_algae_R <- data.frame(sponge_mrm_algae[2])
sponge_mrm_algae_F <- data.frame(sponge_mrm_algae[3])



### Gorgonian ###
gorg_mrm_best <- MRM(gorgonian.bray ~ depth_distance + exposure_dist + linear_dist + mg7_dist, nperm=9999)
gorg_mrm_coral <- MRM(gorgonian.bray ~ depth_distance + exposure_dist + linear_dist + mg7_dist + stony.coral.bray, nperm=9999)
gorg_mrm_algae <- MRM(gorgonian.bray ~ depth_distance + exposure_dist + linear_dist + mg7_dist + macroalgae.bray, nperm=9999)
gorg_mrm_sponge <- MRM(gorgonian.bray ~ depth_distance + exposure_dist + linear_dist + mg7_dist + sponge.bray, nperm=9999)

gorg_mrm_best # 0.62
gorg_mrm_coral # 0.63, coral is significant
gorg_mrm_algae # 0.65, algae is significant
gorg_mrm_sponge # sponge is not significant

# Pull out the sig ones for rmarkdown results
gorg_mrm_coral_coefs <-data.frame(gorg_mrm_coral[1])[6,]
gorg_mrm_coral_R <- data.frame(gorg_mrm_coral[2])
gorg_mrm_coral_F <- data.frame(gorg_mrm_coral[3])

gorg_mrm_algae_coefs <-data.frame(gorg_mrm_algae[1])[6,]
gorg_mrm_algae_R <- data.frame(gorg_mrm_algae[2])
gorg_mrm_algae_F <- data.frame(gorg_mrm_algae[3])



# Save results for rmarkdown read in
corr.mrm.coefs <- list(coral_mrm_algae_coefs,coral_mrm_sponge_coefs,algae_mrm_coral_coefs,
	algae_mrm_sponge_coefs,algae_mrm_gorg_coefs, sponge_mrm_coral_coefs,sponge_mrm_algae_coefs,
	gorg_mrm_coral_coefs,gorg_mrm_algae_coefs)
corr.mrm.R <- list(coral_mrm_algae_R,coral_mrm_sponge_R,algae_mrm_coral_R,
	algae_mrm_sponge_R,algae_mrm_gorg_R, sponge_mrm_coral_R,sponge_mrm_algae_R,
	gorg_mrm_coral_R,gorg_mrm_algae_R)
corr.mrm.F <- list(coral_mrm_algae_F,coral_mrm_sponge_F,algae_mrm_coral_F,
	algae_mrm_sponge_F,algae_mrm_gorg_F, sponge_mrm_coral_F,sponge_mrm_algae_F,
	gorg_mrm_coral_F,gorg_mrm_algae_F)

save(corr.mrm.coefs , file="results/beta-diversity/correlation-mrm-environment-results-coefs.Rdata")
save(corr.mrm.R, file="results/beta-diversity/correlation-mrm-environment-results-R.Rdata")
save(corr.mrm.F, file="results/beta-diversity/correlation-mrm-environment-results-F.Rdata")

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
group <- c("coral", "coral", "coral", "coral", "coral")
coral.plot$group <- group

# Algae
algae.plot = coefDATA(algae_mrm_best)
group <- c("algae", "algae", "algae", "algae")
algae.plot$group <- group

# Sponge
sponge.plot = coefDATA(sponge_mrm_best)
group <- c("sponge", "sponge")
sponge.plot$group <- group

# Gorgonian
gorg.plot = coefDATA(gorg_mrm_best)
group <- c("gorgonian", "gorgonian","gorgonian","gorgonian","gorgonian")
gorg.plot$group <- group

coefs.env <- rbind(coral.plot, algae.plot, sponge.plot, gorg.plot)
coefs.env

# add column for label
coral.plot$label <- "a)"
algae.plot$label <- "b)"
sponge.plot$label <- "c)"
gorg.plot$label <- "d)"

coefs.env <- rbind(coral.plot, algae.plot, sponge.plot, gorg.plot)
head(coefs.env)

coefs.env <- coefs.env %>% filter(env != "Int")


# change order for the plot
coefs.env$group = factor(coefs.env$group, levels=c('coral','algae','sponge','gorgonian'))

###############################################################


###############################################################
# PLOT! 

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
	#legend.position = c(.98,.98),
	#legend.justification= c("right", "top"),
	plot.title=element_text(face= "bold")) + 
  labs(colour = "Environmental Beta Diversity", 
       x= "Environmental Beta Diversity", 
       y= "Beta Diversity Regression Coefficients",
       subtitle = "Abundance") +
	scale_colour_discrete(breaks=c("depth_distance", "exposure_dist", "linear_dist", "mg7_dist", "productivity_dist", "scaridae_dist"),
		labels=c("Depth", "Exposure", "Linear Distance", "Market Gravity", "Productivity", "Grazing")) + 
	coord_flip() + geom_hline(yintercept = 0, linetype=2)+
  ggtitle("Environment")
p

my_tag <- c("(e) Coral", "(f) Algal", "(g) Sponge", "(h) Gorgonian")
p = tag_facet(p, tag_pool = my_tag,
              open = "", close = "",
              fontface = 1,
              hjust = -0.15,x=Inf, y=-Inf)
p

ggsave(file='figures/montastraea/beta-diversity/mrm-env-coefs.pdf', 
       plot=p,
       height=12, width=8)

mrm_bray<- p
save(mrm_bray, file='data/figures/mrm_bray_env.Rdata')
###############################################################




###############################################################
# Pull out the coefficients and p-values of the best env model to plot
# only want to plot the significant ones


# Coral
# algae and sponge are significant
coral.algae.plot = coefDATA(coral_mrm_algae) # function
coral.algae.plot <- coral.algae.plot[nrow(coral.algae.plot), ] # pull out only the last row 
group <- c("coral")
coral.algae.plot$group <- group

coral.sponge.plot = coefDATA(coral_mrm_sponge) # function
coral.sponge.plot <- coral.sponge.plot[nrow(coral.sponge.plot), ] # pull out only the last row 
group <- c("coral")
coral.sponge.plot$group <- group


# Algae
# all three are significant
algae.coral.plot = coefDATA(algae_mrm_coral) # function
algae.coral.plot <- algae.coral.plot[nrow(algae.coral.plot), ] # pull out only the last row 
group <- c("algae")
algae.coral.plot$group <- group

algae.sponge.plot = coefDATA(algae_mrm_sponge) # function
algae.sponge.plot <- algae.sponge.plot[nrow(algae.sponge.plot), ] # pull out only the last row 
group <- c("algae")
algae.sponge.plot$group <- group

algae.gorg.plot = coefDATA(algae_mrm_gorg) # function
algae.gorg.plot <- algae.gorg.plot[nrow(algae.gorg.plot), ] # pull out only the last row 
group <- c("algae")
algae.gorg.plot$group <- group

# Sponge
# coral and algae are significant 
sponge.coral.plot = coefDATA(sponge_mrm_coral) # function
sponge.coral.plot <- sponge.coral.plot[nrow(sponge.coral.plot), ] # pull out only the last row 
group <- c("sponge")
sponge.coral.plot$group <- group

sponge.algae.plot = coefDATA(sponge_mrm_algae) # function
sponge.algae.plot <- sponge.algae.plot[nrow(sponge.algae.plot), ] # pull out only the last row 
group <- c("sponge")
sponge.algae.plot$group <- group

# Gorgonian
# coral and algae are significant
gorg.coral.plot = coefDATA(gorg_mrm_coral) # function
gorg.coral.plot <- gorg.coral.plot[nrow(gorg.coral.plot), ] # pull out only the last row 
group <- c("gorgonian")
gorg.coral.plot$group <- group

gorg.algae.plot = coefDATA(gorg_mrm_algae) # function
gorg.algae.plot <- gorg.algae.plot[nrow(gorg.algae.plot), ] # pull out only the last row 
group <- c("gorgonian")
gorg.algae.plot$group <- group

coefs.taxa <- rbind(coral.algae.plot, coral.sponge.plot, algae.coral.plot, algae.sponge.plot,
	algae.gorg.plot, sponge.coral.plot,sponge.algae.plot, gorg.coral.plot, gorg.algae.plot )
coefs.taxa

# change order for the plot
coefs.taxa$group = factor(coefs.taxa$group, levels=c('coral','algae','sponge','gorgonian'))

# Add the labeler 
coefs.taxa <-coefs.taxa %>% mutate(label=ifelse(group== "coral", "a)", 
                                   ifelse(group== "algae", "b)", 
                                          ifelse(group=='sponge', "c)", "d)"))))
###############################################################




###############################################################
# PLOT! 
p = ggplot(coefs.taxa, aes(x=env, y=coef, color=env)) + geom_point(size=8, shape=18) + facet_grid(group ~ .) + 
  theme(panel.grid.major = element_blank(), 
        panel.border = element_rect(colour = "black", fill= NA), 
        panel.grid.minor = element_blank(), 
		panel.background= element_blank(), 
		axis.text=element_text(size=11), 
		axis.title=element_text(size=13), 
	 	axis.text.y=element_blank(),
		strip.text = element_text(size = 11),
		legend.position= "bottom",
		legend.box= "vertical",
		#legend.text= element_text(size=12),
		#legend.title= element_text(size=12),
		#legend.position = c(.98,.98),
		#legend.justification= c("right", "top"),
		plot.title=element_text(face= "bold")) + 
  labs(colour = "Taxon Beta Diversity", 
       x= "Taxon Beta Diversity",
       y= "  \n   ",
       subtitle= "Abundance") +
	scale_colour_manual(values=c("purple3", "darkgreen", "darkorange2","navyblue"),
		breaks=c("stony.coral.bray","macroalgae.bray", "sponge.bray", "gorgonian.bray"),
		labels=c("Coral", "Algal", "Sponge", "Gorgonian")) + 
	coord_flip() + 
  geom_hline(yintercept = 0, linetype=2)+
  ggtitle("Taxa") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

p

my_tag <- c("(a) Coral", "(b) Algal", "(c) Sponge", "(d) Gorgonian")
p = tag_facet(p, tag_pool = my_tag,
              open = "", close = "",
              fontface= 1,
              hjust = c(-0.65, -0.65, -0.5, -0.39),x=Inf, y=-Inf)
p


mrm_taxa <- p
save(mrm_taxa, file='data/figures/mrm_taxa.Rdata')

ggsave(file='figures/montastraea/beta-diversity/mrm-taxa-coefs.pdf', 
       plot=p,
       height=8, width=7)
###############################################################
