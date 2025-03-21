##mediation analysis of MH EPDS and HAMD with biomarkers
## BLUE = = 500 mg, PINK 1000 mg
##MH outcomes are change scores from baseline, treat as continuous with normal dist

#---- get_data ----

rm(list=ls())
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


##mediator
#results are 500mg - 1000mg at time 2 i.e. 12 months

#### ham_dep

#---- hamd_dep ----
#  25(oh)d3

hamd_dep_25ohd3 = mh_mediate("mh_hamd_dep", "bl_mh_hamd_dep", 
                   "x25ohd3fu", "x25ohd3_screen",  nsim = 500)
summary(hamd_dep_25ohd3)

plot(hamd_dep_25ohd3, treatment = "both")




# 25(oh)d3 at 6 weeks only

hamd_dep_25ohd3_fu1 = mh_mediate_lm("mh_hamd_dep", "bl_mh_hamd_dep", 
                              "x25ohd3fu", "x25ohd3_screen", timept = 1, nsim = 500)
summary(hamd_dep_25ohd3_fu1)

plot(hamd_dep_25ohd3_fu1, treatment = "both")

#sensitivity, only available for one time point
sens.cont <- medsens(hamd_dep_25ohd3_fu1, rho.by=.1, eps=.01, effect.type="both")

# Use summary function to display results
summary(sens.cont)

# Plot true ACMEs and ADEs as functions of rho
par.orig <- par(mfrow = c(2,1))
plot(sens.cont,  ylim=c(-.2,.2))


## il-6
hamd_dep_il6 = mh_mediate ("mh_hamd_dep", "bl_mh_hamd_dep", 
                           "log_il6", "log_il6_screen", nsim = 500)

summary(hamd_dep_il6)
plot(hamd_dep_il6)


## il-6, 6 weeks only
hamd_dep_il6_fu1 = mh_mediate_lm("mh_hamd_dep", "bl_mh_hamd_dep", 
                           "log_il6", "log_il6_screen", timept = 1, nsim = 500)

summary(hamd_dep_il6_fu1)
plot(hamd_dep_il6_fu1)

sens.cont_il6 <- medsens(hamd_dep_il6_fu1, rho.by=.1, eps=.01, effect.type="both")

# Use summary function to display results
summary(sens.cont_il6)

# Plot true ACMEs and ADEs as functions of rho
plot(sens.cont_il6,  ylim=c(-.2,.2))

par(mfrow = c(1,1))

#---- epds ----
# 25(oh)d3

epds_25ohd3 = mh_mediate ("mh_epds", "bl_mh_epds", 
                           "x25ohd3fu", "x25ohd3_screen", nsim = 500)
summary(epds_25ohd3)
plot(epds_25ohd3)


# 6 weeks only

epds_25ohd3_fu1 = mh_mediate_lm("mh_epds", "bl_mh_epds", 
                                    "x25ohd3fu", "x25ohd3_screen", timept = 1, nsim = 500)
summary(epds_25ohd3_fu1)

plot(epds_25ohd3_fu1)


## il-6

epds_il6 = mh_mediate ("mh_epds", "bl_mh_epds", 
                       "log_il6", "log_il6_screen", nsim = 500)

summary(epds_il6)
plot(epds_il6)






#---- lmer_mh_cumdose ----
#### MH outcomes

mh <- mh %>% mutate(cumdose = cumdose/500)

MH_cumdose <- function(outvar, exposure = "cumdose", 
                       adjvars = c("Age","income","Education_years"),
                       basevar = NULL, exp = F) {
  
  form = formula(paste0(outvar, " ~",  exposure,"  + " , paste(adjvars, collapse = "+"),
                        "+ (1|Study_ID) +",basevar) )
  
  fit <- lmer(formula = form,data = mh)
  s <- summary(fit)$coefficients 
  
  df <- data.frame(Outcome = outvar,
                   Estimate = s[exposure,"Estimate"],
                   se = s[exposure,"Std. Error"],
                   p = format.pval(s[exposure,"Pr(>|t|)"], eps=0.001, 2)) %>%
    mutate(lower = Estimate - 1.96*se,
           upper = Estimate + 1.96*se)
  
  if (exp) {
    
    df <- df %>%
      mutate(est_ci = paste0(round(exp(Estimate), 2), " (",
                             round(exp(lower), 2), ", ",
                             round(exp(upper), 2), ")"))
  } else {
    
    df <- df %>%
      mutate(
        est_ci = paste0(round(Estimate, 2), " (",
                        round(lower, 2), ", ",
                        round(upper, 2), ")"))
    
  }
  
  
  return(df)
  
}

#cumulative dose

(hamd <- MH_cumdose(outvar = "MH_HAMD_DEP" ,  basevar = "BL_MH_HAMD_DEP"))

(hama <- MH_cumdose(outvar = "MH_HAMA_ANX",  
                    basevar = "BL_MH_HAMA_ANX")) 

(SIGHAD <- MH_cumdose(outvar = "MH_SIGHAD", basevar ="BL_MH_SIGHAD"))


(GAF <- MH_cumdose(outvar = "MH_GAF", basevar ="BL_MH_GAF"))

(FAST <- MH_cumdose(outvar = "MH_FAST", basevar ="BL_MH_FAST"))

(epds <- MH_cumdose(outvar = "MH_EPDS", basevar ="BL_MH_EPDS"))



cumdose_tab <- bind_rows(hamd, hama, SIGHAD, GAF, FAST, epds) %>%
  mutate(Outcome = gsub("MH_","",Outcome)) %>%
  select(Outcome, est_ci,p)


cumdose_tab %>%
  flextable %>%
  autofit() %>%
  theme_box() %>%
  set_header_labels(est_ci = "Estimate (increase in 500mg) (95% CI)")


##topups
(hamd_top <- MH_cumdose(outvar = "MH_HAMD_DEP",  basevar = "BL_MH_HAMD_DEP", 
                        adjvars = c("Age","income","Education_years","rand"),
                        exposure = "ntopup", exp = F))

(hama_top <- MH_cumdose(outvar = "MH_HAMA_ANX", basevar = "BL_MH_HAMA_ANX", 
                        adjvars = c("Age","income","Education_years","rand"),
                        exposure = "ntopup", exp = F)) 


(SIGHAD_top <- MH_cumdose(outvar = "MH_SIGHAD", basevar ="BL_MH_SIGHAD", 
                          adjvars = c("Age","income","Education_years","rand"),
                          exposure = "ntopup", exp = F))


(GAF_top <- MH_cumdose(outvar = "MH_GAF",basevar = "BL_MH_GAF", 
                       adjvars = c("Age","income","Education_years","rand"),
                       exposure = "ntopup", exp = F))

(FAST_top <- MH_cumdose(outvar = "MH_FAST", basevar ="BL_MH_FAST", 
                        adjvars = c("Age","income","Education_years","rand"),
                        exposure = "ntopup", exp = F))

(epds_top <- MH_cumdose("MH_EPDS", basevar ="BL_MH_EPDS", 
                        adjvars = c("Age","income","Education_years","rand"),
                        exposure = "ntopup", exp = F))


ntopup <- bind_rows(hamd_top, hama_top, SIGHAD_top, GAF_top, FAST_top, epds_top) %>%
  mutate(Outcome = gsub("MH_","",Outcome)) %>%
  select(Outcome, est_ci,p)


ntopup %>%
  flextable %>%
  autofit() %>%
  theme_box() %>%
  set_header_labels(est_ci = "Estimate (every extra topup) (95% CI)")


##add topup to secondary outcome models

MH_topup <- function(outvar, basevar) {
  
  form = formula(paste0(outvar, " ~ rand*time + stratification + Age + income + 
                        Education_years + ntopup + (1|Study_ID)")) #
  
  
  fit <- lmer(formula = form,data = mh)
  
  
  rt <- emmeans(fit, ~ rand * time )  #estimated means
  (pairs <- summary(pairs(rt, simple="rand")))  #contrasts
  
  
  estimate3<-as.numeric(pairs[1,"estimate"])  #estimate  at 6 weeks
  se3<-as.numeric(pairs[1,"SE"]) #se at 6 weeks
  p3 <-pairs[1,"p.value"]
  
  estimate6<-as.numeric(pairs[2,"estimate"])  #estimate at 12m
  se6<-as.numeric(pairs[2,"SE"]) #se at 12m
  p6 <-pairs[2,"p.value"]
  
  
  
  Outcome=c(outvar,"")
  time=c("6 weeks", "12 months")
  estimate = rbind(estimate3,estimate6)
  se=rbind(se3,se6)
  p=rbind(p3,p6)
  
  df = data.frame(Outcome=Outcome, time=time, Estimate = estimate, 
                  se = se, p=p)
  
  rownames(df) <- c()
  return(df)
}


(hamd_adj <- MH_topup(outvar = "MH_HAMD_DEP" ,  basevar = "BL_MH_HAMD_DEP"))

(hama_adj <- MH_topup(outvar = "MH_HAMA_ANX",  
                      basevar = "BL_MH_HAMA_ANX")) 

(SIGHAD_adj <- MH_topup(outvar = "MH_SIGHAD", basevar ="BL_MH_SIGHAD"))


(GAF_adj <- MH_topup(outvar = "MH_GAF", basevar ="BL_MH_GAF"))

(FAST_adj <- MH_topup(outvar = "MH_FAST", basevar ="BL_MH_FAST"))

(epds_adj <- MH_topup("MH_EPDS", basevar ="BL_MH_EPDS"))

mh_adj <- bind_rows(hamd_adj, hama_adj, SIGHAD_adj, GAF_adj, FAST_adj, epds_adj) %>%
  mutate(Outcome = gsub("MH_","",Outcome),
         est_ci = paste0(round(Estimate, 2), " (",
                         round(Estimate - 1.96*se, 2), ", ",
                         round(Estimate + 1.96*se, 2), ")"),
         p=format.pval(p, eps=0.001, 2)) %>%
  select(Outcome, time, est_ci,p)

mh_adj %>%
  flextable %>%
  autofit() %>%
  theme_box() %>%
  set_header_labels(est_ci = "500mg - 1000mg (95% CI)",
                    time= "Time Point")


