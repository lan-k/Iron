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
  mutate(treat = factor(ifelse(rand == "BLUE", 1 ,0 )), #500mg is treat = 1
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

# 25(oh)d3 at 6 weeks
hamd_dep_25ohd3_fu1 = mh_mediate(outcome = "mh_hamd_dep", outcome_base ="bl_mh_hamd_dep", 
                             mediator= "x25ohd3fu", mediator_base = "x25ohd3_screen", 
                             timept=1,
                             nsim = 500)
df_hamd_dep_25ohd3_fu1= tab_mh_mediate(hamd_dep_25ohd3_fu1$medfit, "HAM-D", "25(OH)D3")

# plot(hamd_dep_25ohd3_fu1$medfit, treatment = "both", main = "HAMD and 25(OH)D3 at 6 weeks",
#      xlab = "Difference in HAMD change from baseline (500mg-1000mg)")

# 25(oh)d3 at 12 months
hamd_dep_25ohd3_fu2 = mh_mediate(outcome = "mh_hamd_dep", outcome_base ="bl_mh_hamd_dep", 
                                 mediator= "x25ohd3fu", mediator_base = "x25ohd3_screen", 
                                 timept=2,
                                 nsim = 500)
df_hamd_dep_25ohd3_fu2= tab_mh_mediate(hamd_dep_25ohd3_fu2$medfit, "HAM-D", "25(OH)D3")

# plot(hamd_dep_25ohd3_fu2$medfit, treatment = "both", main = "HAMD and 25(OH)D3 at 12 months",
#      xlab = "Difference in HAMD change from baseline (500mg-1000mg)")


df_hamd_dep_25ohd3 = bind_rows(df_hamd_dep_25ohd3_fu1 %>% mutate(timept = "6 weeks"),
                               df_hamd_dep_25ohd3_fu2 %>% 
                                 mutate(timept = "12 months",
                                        Outcome="", Mediator = "") )

# 25(oh)d3 at 6 weeks only
# hamd_dep_25ohd3_fu1 = mh_mediate_lm("mh_hamd_dep", "bl_mh_hamd_dep", 
#                               "x25ohd3fu", "x25ohd3_screen", timept = 1, nsim = 500)
# 
# 
# df_hamd_dep_25ohd3_fu1= tab_mh_mediate(hamd_dep_25ohd3_fu1, 
#                                    "HAMD_DEP at 6 weeks", "25(OH)D3")
# 
# 
# plot(hamd_dep_25ohd3_fu1, treatment = "both", main = "HAMD and 25(OH)D3 at 6 weeks",
#      xlab = "Difference in HAMD change from baseline (500mg-1000mg)")

# epds

# 6 weeks

epds_25ohd3_fu1 = mh_mediate ("mh_epds", "bl_mh_epds", 
                          "x25ohd3fu", "x25ohd3_screen",timept=1, nsim = 500)
df_epds_25ohd3_fu1= tab_mh_mediate(epds_25ohd3_fu1$medfit, "EPDS", "25(OH)D3")

# plot(epds_25ohd3_fu1$medfit,treatment = "both", main = "EPDS and 25(OH)D3 at 6 weeks",
#      xlab = "Difference in EPDS change from baseline (500mg-1000mg)",
#      xlim=c(-4,4))



# 12 months

epds_25ohd3_fu2 = mh_mediate("mh_epds", "bl_mh_epds", 
                            "x25ohd3fu", "x25ohd3_screen", timept = 2, nsim = 500)


df_epds_25ohd3_fu2= tab_mh_mediate(epds_25ohd3_fu2$medfit, "EPDS", "25(OH)D3")
# plot(epds_25ohd3_fu2$medfit, treatment = "both", main = "EPDS and 25(OH)D3 at 12 months",
#      xlab = "Difference in EPDS change from baseline (500mg-1000mg)",
#      xlim=c(-4,4))



df_epds_25ohd3 = bind_rows(df_epds_25ohd3_fu1 %>% mutate(timept = "6 weeks"),
                               df_epds_25ohd3_fu2 %>% 
                             mutate(timept = "12 months",
                                    Outcome="", Mediator = "") )

#---- il6 ----

## hamd

# 6 weeks

hamd_dep_il6_fu1 = mh_mediate("mh_hamd_dep", "bl_mh_hamd_dep","log_il6", "log_il6_screen", 
                                  timept = 1, exp=T, nsim = 500)

df_hamd_dep_il6_fu1= tab_mh_mediate(hamd_dep_il6_fu1$medfit, "HAM-D", "IL-6")
# plot(hamd_dep_il6_fu1$medfit, treatment = "both", main = "HAMD and IL-6 at 6 weeks",
#      xlab = "Difference in HAMD change from baseline (500mg-1000mg)", xlim=c(-5,1))
# 

#12 months
hamd_dep_il6_fu2 = mh_mediate("mh_hamd_dep", "bl_mh_hamd_dep", 
                              "log_il6", "log_il6_screen", timept = 2,exp=T, nsim = 500)

df_hamd_dep_il6_fu2= tab_mh_mediate(hamd_dep_il6_fu2$medfit, "HAM-D", "IL-6")


# plot(hamd_dep_il6_fu2$medfit, treatment = "both", main = "HAMD and IL-6",
#      xlab = "Difference in HAMD change from baseline (500mg-1000mg)", xlim=c(-5,1))


df_hamd_dep_il6 = bind_rows(df_hamd_dep_il6_fu1 %>% mutate(timept = "6 weeks"),
                               df_hamd_dep_il6_fu2 %>% 
                              mutate(timept = "12 months",
                                     Outcome="", Mediator = ""))


## epds
# 6 weeks
epds_il6_fu1 = mh_mediate("mh_epds", "bl_mh_epds", 
                       "log_il6", "log_il6_screen", timept=1, exp=T, nsim = 500)


df_epds_il6_fu1= tab_mh_mediate(epds_il6_fu1$medfit, "EPDS", "IL-6")


# plot(epds_il6_fu1$medfit, treatment = "both", main = "EPDS and IL-6 at 6 weeks",
#      xlab = "Difference in EPDS change from baseline (500mg-1000mg)")
# 



# 12 months
epds_il6_fu2 = mh_mediate("mh_epds", "bl_mh_epds", 
                             "log_il6", "log_il6_screen", timept = 2, exp=T, nsim = 500)


df_epds_il6_fu2= tab_mh_mediate(epds_il6_fu2$medfit, "EPDS", "IL-6")


# plot(epds_il6_fu2$medfit,treatment = "both", main = "EPDS and IL-6 at 12 months",
#      xlab = "Difference in EPDS change from baseline (500mg-1000mg)")


df_epds_il6 = bind_rows(df_epds_il6_fu1 %>% mutate(timept = "6 weeks"),
                        df_epds_il6_fu2 %>% 
                          mutate(timept = "12 months",
                                 Outcome="", Mediator = "") )
