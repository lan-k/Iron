##mediation analysis of MH EPDS and HAMD with biomarkers
## BLUE = = 500 mg, PINK 1000 mg
##MH outcomes are change scores from baseline, treat as continuous with normal dist

#results are 500mg - 1000mg at time 2 i.e. 12 months

#---- get_data ----

rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(mediation)
  library(lmerTest)
  library(lme4)
  library(emmeans)
  library(splines)
  library(flextable)
  library(visdat)
})



load(file="mh_il6_vitd.rds" )
source("99_functions.R")


mh_il6_vitd <- mh_il6_vitd %>%
  mutate(treat = factor(ifelse(rand == "BLUE", 1 ,0 )),
         il_6_nolim = ifelse(il_6 == 0.2, 0.2*runif(1),il_6),
         il_6_nolim_screen = ifelse(il_6_screen == 0.2, 0.2*runif(1),il_6_screen),
         log_il6 = log(il_6_nolim),
         log_il6_screen = log(il_6_nolim_screen))

##check missing data pattern

mh_il6_vitd %>%
  # group_by(time) %>%
  dplyr::select(study_id, income,
                education_years, x25ohd3fu, x25ohd3_screen, il_6, il_6_screen,
                mh_hamd_dep, bl_mh_hamd_dep,
                mh_epds, bl_mh_epds)  %>% 
  visdat::vis_miss()




#---- 25ohd3 ----
##mediator 25(oh)d3
#### ham_dep


hamd_dep_25ohd3 = mh_mediate("mh_hamd_dep", "bl_mh_hamd_dep", 
                   "x25ohd3fu", "x25ohd3_screen",  nsim = 500)
df_hamd_dep_25ohd3= tab_mh_mediate(hamd_dep_25ohd3$medfit, "HAMD_DEP", "25(OH)D3")

plot(hamd_dep_25ohd3$medfit, treatment = "both", main = "HAMD and 25(OH)D3 over 12 months",
     xlab = "Difference in HAMD change from baseline (500mg-1000mg)")




# 25(oh)d3 at 6 weeks only

hamd_dep_25ohd3_fu1 = mh_mediate_lm("mh_hamd_dep", "bl_mh_hamd_dep", 
                              "x25ohd3fu", "x25ohd3_screen", timept = 1, nsim = 500)


df_hamd_dep_25ohd3_fu1= tab_mh_mediate(hamd_dep_25ohd3_fu1, 
                                   "HAMD_DEP at 6 weeks", "25(OH)D3")


plot(hamd_dep_25ohd3_fu1, treatment = "both", main = "HAMD and 25(OH)D3 at 6 weeks",
     xlab = "Difference in HAMD change from baseline (500mg-1000mg)")

# #sensitivity, only available for one time point
# sens.cont <- medsens(hamd_dep_25ohd3_fu1, rho.by=.1, eps=.01, effect.type="both")
# 
# # Use summary function to display results
# summary(sens.cont)
# 
# # Plot true ACMEs and ADEs as functions of rho
# par.orig <- par(mfrow = c(2,1))
# plot(sens.cont,  ylim=c(-.2,.2))

# epds

epds_25ohd3 = mh_mediate ("mh_epds", "bl_mh_epds", 
                          "x25ohd3fu", "x25ohd3_screen", nsim = 500)
df_epds_25ohd3= tab_mh_mediate(epds_25ohd3$medfit, "EPDS", "25(OH)D3")

plot(epds_25ohd3$medfit, treatment = "both", main = "EPDS and 25(OH)D3 over 12 months",
     xlab = "Difference in EPDS change from baseline (500mg-1000mg)")

# 6 weeks only

epds_25ohd3_fu1 = mh_mediate_lm("mh_epds", "bl_mh_epds", 
                                "x25ohd3fu", "x25ohd3_screen", timept = 1, nsim = 500)


df_epds_25ohd3_fu1= tab_mh_mediate(epds_25ohd3_fu1, "EPDS", "25(OH)D3")


plot(epds_25ohd3_fu1,treatment = "both", main = "EPDS and 25(OH)D3 at 6 weeks",
     xlab = "Difference in EPDS change from baseline (500mg-1000mg)")



desc_25ohd3-> 


#---- il6 ----

## hamd
hamd_dep_il6 = mh_mediate ("mh_hamd_dep", "bl_mh_hamd_dep", 
                           "log_il6", "log_il6_screen", nsim = 500)

df_hamd_dep_il6= tab_mh_mediate(hamd_dep_il6$medfit, "HAMD_DEP at 6 weeks", "IL-6")


plot(hamd_dep_il6$medfit, treatment = "both", main = "HAMD and IL-6 over 12 months",
     xlab = "Difference in HAMD change from baseline (500mg-1000mg)")


## il-6, 6 weeks only
hamd_dep_il6_fu1 = mh_mediate_lm("mh_hamd_dep", "bl_mh_hamd_dep", 
                                 "log_il6", "log_il6_screen", timept = 1, nsim = 500)

summary(hamd_dep_il6_fu1)
plot(hamd_dep_il6_fu1)

# sens.cont_il6 <- medsens(hamd_dep_il6_fu1, rho.by=.1, eps=.01, effect.type="both")
# 
# # Use summary function to display results
# summary(sens.cont_il6)
# 
# # Plot true ACMEs and ADEs as functions of rho
# plot(sens.cont_il6,  ylim=c(-.2,.2))
# 
# par(mfrow = c(1,1))


## epds

epds_il6 = mh_mediate("mh_epds", "bl_mh_epds", 
                       "log_il6", "log_il6_screen", nsim = 500)


df_epds_il6= tab_mh_mediate(epds_il6$medfit, "EPDS", "IL-6")

plot(epds_il6$medfit, treatment = "both", main = "EPDS and IL-6 over 12 months",
     xlab = "Difference in EPDS change from baseline (500mg-1000mg)")

plot(epds_il6)


# 6 weeks only
epds_il6_fu1 = mh_mediate_lm("mh_epds", "bl_mh_epds", 
                             "log_il6", "log_il6_screen", timept = 1, nsim = 500)


df_epds_il6_fu1= tab_mh_mediate(epds_il6_fu1, "EPDS", "IL-6")


plot(epds_il6_fu1,treatment = "both", main = "EPDS and IL-6 at 6 weeks",
     xlab = "Difference in EPDS change from baseline (500mg-1000mg)")


