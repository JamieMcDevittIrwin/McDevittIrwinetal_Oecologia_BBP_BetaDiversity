### Put both environment and taxa beta-div plots together ###
############################################################

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
library(gridExtra)
library(ggpubr)

load('data/figures/mrm_bray_env.Rdata')
load('data/figures/mrm_sorenson_env.Rdata')
load('data/figures/mrm_taxa.Rdata')
ls()


# Arrange the two environmental beta diversity plots together
p <- ggpubr::ggarrange(mrm_bray, mrm_sorenson, 
          common.legend=TRUE, 
          legend='bottom')
p

ggsave(p, file="figures/montastraea/beta-diversity/combined-beta-env.pdf",
       width=10, height=8)


# Combine environmental beta diversity with taxa beta diversity
p_mrmtaxa <- mrm_taxa

p_full <- grid.arrange(p_mrmtaxa, 
             p, 
             ncol=3, nrow=1,
             layout_matrix= rbind(c(1,2,2)))
p_full

ggsave(p_full, file="figures/montastraea/beta-diversity/combined-beta-env-taxa.pdf",
       width=16, height=12)
 


