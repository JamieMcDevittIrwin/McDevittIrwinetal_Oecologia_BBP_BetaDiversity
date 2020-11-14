# Plot Best Models for Species richness

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
library(ggplot2)
library(sandwich)
library(lmtest)
library(egg)

# load in data
richness.metadata.sc <- read.csv("data/calculated/richness.metadata.scaled.csv")


# UPDATE: sponge best model now has glm instead of lm, using robust SE's from sandwich package for all models to be conservative on heteroskedasticity

# Best Model (determined from scripts)
# Coral: lm(richness_coral ~ Net_primary_productivity + Exposure, data = richness.metadata.sc)
# Algae: lm(richness_algae ~ Depth, data = richness.metadata.sc)
# Sponge: glm(richness_sponge ~ Depth + MG...Bahamas.big.7 + Net_primary_productivity,family=poisson, data=richness.metadata.sc)
# Gorgonian: lm(richness_gorgonian ~ Net_primary_productivity + Exposure + Scaridae + MG...Bahamas.big.7,  data = richness.metadata.sc)
#######################################################################



#######################################################################
# EXTRACT FOR PLOT # 

# Best models
M.coral <- lm(richness_coral ~ Net_primary_productivity + Exposure, data = richness.metadata.sc)
M.sponge <- glm(richness_sponge ~ Depth + MG...Bahamas.big.7 + Net_primary_productivity,family=poisson, data=richness.metadata.sc) 
M.gorg <- lm(richness_gorgonian ~ Net_primary_productivity 
	+ Exposure + Scaridae + MG...Bahamas.big.7,  data = richness.metadata.sc)

# add name of variable so the function knows which explan var to pull out
coral.fix <- c("Net_primary_productivity", "Exposure")
sponge.fix <- c("Depth", "MG...Bahamas.big.7", "Net_primary_productivity")
gorg.fix <- c("Net_primary_productivity", "Exposure", "Scaridae", "MG...Bahamas.big.7")



# Extra slopes and SEs #
# test how to extract first
coeftest(M.gorg, vcov= vcovHC)["Net_primary_productivity", 1] # pulls out coef
coeftest(M.gorg, vcov= vcovHC)["Net_primary_productivity", 2] # pulls out SE

# now write a function based on this
plotgg <- function(MODEL, FIX){
	#extract slopes and SE
	slopes <- coeftest(MODEL, vcov= vcovHC)[FIX, 1]  #' slopes
	ses <- coeftest(MODEL, vcov= vcovHC)[FIX, 2]  #' SEs

	# convert to df
	data= data.frame(ses,slopes)
	data <- tibble::rownames_to_column(data)
	return(data)
}


coral.df <- plotgg(M.coral, coral.fix)
sponge.df <- plotgg(M.sponge, sponge.fix)
gorg.df <- plotgg(M.gorg, gorg.fix)



# add column for each group to facet 
coral.df$facet <- "coral"
sponge.df$facet <- "sponge"
gorg.df$facet <- "gorgonian"

full.df <- rbind(coral.df, sponge.df, gorg.df)
head(full.df)

# add column for label
coral.df$label <- "a)"
sponge.df$label <- "b)"
gorg.df$label <- "c)"

full.df <- rbind(coral.df, sponge.df, gorg.df)
full.df

# change order for the plot
full.df$facet = factor(full.df$facet, levels=c('coral','sponge','gorgonian'))
head(full.df)

#######################################################################


#######################################################################
# PLOT # 
p = ggplot(full.df, aes(x=rowname, y=slopes, color=rowname)) + 
  geom_pointrange(aes(ymin=slopes-ses, ymax=slopes+ses), size=1.5)  + 
  facet_grid(facet ~ ., scales= "free_y") + theme(
	panel.grid.major = element_blank(), 
	panel.border = element_rect(colour = "black", fill= NA), 
	panel.grid.minor = element_blank(), 
		panel.background= element_blank(), 
	axis.text=element_text(size=11), 
	axis.title=element_text(size=13), 
	 	#axis.text.y=element_blank(),	
	legend.position = c(.98,.98),
	legend.justification= c("right", "top"),
	axis.text.x=element_text(angle=25, hjust=0.5, vjust=0.7),
	strip.text = element_text(size = 11)) + 
  labs(colour = "Environmental Drivers", 
       x= "Environmental Drivers", 
       y= "Alpha Diversity Coefficients") + 
	#coord_flip() + 
  geom_hline(yintercept = 0, linetype=2) +
	scale_colour_manual(values=c("#B79F00", "#00BFC4", "#619CFF", "#F564E3", "#F8766D"),
		breaks=c("Exposure", 
		         "MG...Bahamas.big.7", 
		         "Net_primary_productivity", 
		         "Scaridae",
		         "Depth"),
		labels=c("Exposure", 
		         "Market Gravity", 
		         "Productivity", 
		         "Grazing",
		         "Depth")) + 
  scale_x_discrete(labels=c("Exposure" = "Exposure", 
                            "MG...Bahamas.big.7" = "Market Gravity",
                            "Net_primary_productivity" = "Productivity",
                            "Scaridae" = "Grazing", "Depth"= "Depth"))

# labels
my_tag <- c("(a) Coral", "(b) Sponge", "(c) Gorgonian")
p = tag_facet(p, tag_pool = my_tag,
              open = "", close = "",
              fontface= 1, 
              hjust = -0.25)
p


ggsave(filename='figures/montastraea/species_richness/env-coefs.pdf',
       plot=p,
       height=11,
       width=8)


# this shows you the default ggplot colours depending on your n explanatory vars
library(scales)
show_col(hue_pal()(6))

#######################################################################



#######################################################################
# Robust Standard Errors#
library("sandwich")

# Gorgonian
coef(M.gorg)
summary(M.gorg)

# default is HC3
coeftest(M.gorg, vcov= vcovHC) # changes the SE's by 0.03-0.1 and makes scaridae not significant anymore
coeftest(M.gorg)

# Sponge
coef(M.sponge) # gives you the same output as the summary result
summary(M.sponge)
coeftest(M.sponge, vcov= vcovHC(M.sponge, type="HC3")) 
#coeftest(M.sponge, vcov= vcovHC(M.sponge,type="HC0")) 


# Coral
summary(M.coral)
coeftest(M.coral, vcov= vcovHC) # one of the standard errors increased and one decreased (making one sig and one not sig anymore)


