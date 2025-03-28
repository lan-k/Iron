##time varying mediation with tvmediation
# https://cran.r-project.org/web/packages/tvmediation/vignettes/Time_Varying_Mediation_Analysis.html

rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(tvmediation)
library(lmerTest)
library(lme4)
library(emmeans)
library(splines)
library(flextable)
library(visdat)

# citation("bmlm")

load(file="mh_il6_vitd.rds" )
source("99_functions.R")


mh_il6_vitd <- mh_il6_vitd %>%
  mutate(treat = factor(ifelse(rand == "BLUE", 1 ,0 )),
         il_6_nolim = ifelse(il_6 == 0.2, 0.2*runif(1),il_6),
         il_6_nolim_screen = ifelse(il_6_screen == 0.2, 0.2*runif(1),il_6_screen),
         log_il6 = log(il_6_nolim),
         log_il6_screen = log(il_6_nolim_screen)) 

#convert to wide format for the tma function
