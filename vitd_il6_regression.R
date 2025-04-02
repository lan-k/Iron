## lmer and tobit regression
# ---- get_data ----

# rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(lmerTest)
library(emmeans)
library(censReg)
library(plm)
library(splines)
library(flextable)
  library(ggplot2)
})


load(file="mh_il6_vitd.rds" )
source("99_functions.R")

mh_il6_vitd %>%
  # group_by(time) %>%
  dplyr::select(study_id, income,
                education_years, 
                mh_hamd_dep, bl_mh_hamd_dep,
                mh_epds, bl_mh_epds)  %>% 
  na.omit() %>% 
  dplyr::select(study_id, income,
                education_years, x25ohd3fu, x25ohd3_screen, il_6, il_6_screen,
                mh_hamd_dep, bl_mh_hamd_dep,
                mh_epds, bl_mh_epds)



# ---- desc ----

# hist(base$x25ohd3) #most normal
# hist(base$x25ohd2)
# hist(base$epi25ohd3) 


base_n = base %>%
  group_by(rand) %>%
  summarise(n=n()) %>%
  ungroup()


blue_lab = paste0("500 mg"," (n=", as.character(base_n %>% filter(rand=="BLUE") %>% pull(n)), ")")
pink_lab = paste0("1000 mg"," (n=", as.character(base_n %>% filter(rand=="PINK") %>% pull(n)), ")")


base_desc <- base %>%
  group_by(rand) %>%
  reframe(across(matches("il_6|oh"), \(x) mean_log(x, digits=2))) %>% #median_iqr
  ungroup() %>%
  mutate(FU="SCR")


base_desc_wide <- base_desc %>%
  pivot_wider(names_from = "rand", 
              values_from = c("il_6","epi25ohd3","x25ohd2", "x25ohd3") )


desc = mh_il6_vitd %>%
  select(!contains("screen")) %>%
  group_by(rand, FU) %>%
  reframe(across(matches("il_6|oh"), \(x) mean_log(x, digits=2))) %>% #median_iqr
  ungroup()

colnames(desc) <-  sub("fu", "", colnames(desc))



desc_wide <- desc %>%
  pivot_wider(names_from = "rand", 
              values_from = c("il_6","epi25ohd3","x25ohd2", "x25ohd3") )

desc_wide <-bind_rows(base_desc_wide, desc_wide) 



# desc_wide %>%
#   select(FU, il_6_BLUE,il_6_PINK,epi25ohd3_BLUE, epi25ohd3_PINK,
#          x25ohd2_BLUE, x25ohd2_PINK, x25ohd3_BLUE, x25ohd3_PINK) %>%
#   mutate(FU = case_when(FU == "SCR" ~ "Screening",
#                         FU == "FU1" ~ "6 weeks PP",
#                         FU == "FU2" ~ "12 months PP")) %>%
#   flextable() %>%
#   autofit() %>%
#   set_header_labels(FU = "Time point",
#                     il_6_BLUE = "IL-6, \n500 mg",
#                     il_6_PINK = "IL-6, \n1000 mg",
#                     epi25ohd3_BLUE = "Epi-25(OH)D3, \n500mg",
#                     epi25ohd3_PINK = "Epi-25(OH)D3, \n1000mg",
#                     x25ohd2_BLUE = "25(OH)D2, \n500mg",
#                     x25ohd2_PINK = "25(OH)D2, \n1000mg",
#                     x25ohd3_BLUE = "25(OH)D3, \n500mg",
#                     x25ohd3_PINK = "25(OH)D3, \n1000mg"
#                     ) %>%
#   set_caption(caption = "Biomarkers (mean (SD))")
# 


# ---- lmer_25ohd3_il_6 ----
# repeated measures with a random numbe below the lower limit
#Epi 25(OH)D3 lower limit is 2.0
# 25(OH)D2 lower limit is 3.0
# IL-6 lower limit is 0.2

# 25ohd3 not affected by detection limits
# hist(mh_il6_vitd$x25ohd3_screen)
# hist(mh_il6_vitd$x25ohd3fu)



#Epi 25(OH)D3 lower limit is 2.0
# 25(OH)D2 lower limit is 3.0
# IL-6 lower limit is 0.2
mh_il6_vitd <- mh_il6_vitd %>% 
  mutate(id=paste(rand, study_id, sep="_"),
         time=factor(time),
         epi25ohd3fu_cens = ifelse(epi25ohd3fu == 1.99, 2, epi25ohd3fu),
         epi25ohd3_screen_cens = ifelse(epi25ohd3fu == 1.99, 2*runif(1), 
                                        epi25ohd3_screen),
         x25ohd2fu_cens= ifelse(x25ohd2fu == 2.99, 3, x25ohd2fu),
         x25ohd2_screen_cens = ifelse(x25ohd2_screen == 2.99, 3*runif(1), x25ohd2_screen),
         il_6_nolim = ifelse(il_6 == 0.2, 0.2*runif(1),il_6),
         il_6_nolim_screen = ifelse(il_6_screen == 0.2, 0.2*runif(1),il_6_screen),
         x25ohd3_change = x25ohd3fu - x25ohd3_screen) 




# summary(mh_il6_vitd$x25ohd2fu)
# hist(mh_il6_vitd$x25ohd2fu)

# repeated measures for 25(OH)D3

# form = formula(log(x25ohd3fu) ~ rand*time + stratification + log(x25ohd3_screen) +
#                  (1|study_id))
# 
# 
# fit_25ohd3 <- lmer(formula = form, data = mh_il6_vitd)
# summary(fit_25ohd3)
# 
# 
# # anova(fit_25ohd3)
# (p <- summary(pairs(emmeans(fit_25ohd3, ~ rand * time ),
#                     simple="rand")))  #contrasts PINK - BLUE i.e. 1000 - 500
# lsmeans(fit_25ohd3, ~ rand * time)
# 
# estimate3<-as.numeric(p[1,"estimate"])  #estimate  at 6 weeks
# se3<-as.numeric(p[1,"SE"]) #se at 6 weeks
# p3 <-format.pval(p[1,"p.value"], eps=0.001,2)
# 
# estimate6<-as.numeric(p[2,"estimate"])  #estimate at 12m
# se6<-as.numeric(p[2,"SE"]) #se at 12m
# p6 <-format.pval(p[2,"p.value"], eps=0.001,2)
# 
# time=c("6 weeks", "12 months")
# estimate = rbind(estimate3,estimate6)
# se=rbind(se3,se6)
# p=rbind(p3,p6)
# 
# df_25ohd3 = data.frame(Outcome = c("25(OH)D3",""),
#                        Contrast = c("500mg vs 1000mg",""), 
#                        Level=time, Estimate = exp(-estimate), 
#                        lower = exp(-estimate - 1.96*se),
#                        upper = exp(-estimate + 1.96*se),
#                        se = se, p=p)
# 
# rownames(df_25ohd3) <- c()

# change from baseline
# fit_change = lmer(x25ohd3_change ~ rand*time + stratification + log(x25ohd3_screen) +
#                     (1|study_id),data = mh_il6_vitd)
# emmeans(fit_change, ~ rand * time)


# il_6 using a random number between 0 and 0.2 for detection limit


# hist(mh_il6_vitd$il_6_nolim)
# hist(mh_il6_vitd$il_6_nolim_screen)
# hist(log(mh_il6_vitd$il_6_nolim))
# hist(log(mh_il6_vitd$il_6_nolim_screen))


# form = as.formula(log(il_6_nolim) ~ rand*time + stratification + 
#                     log(il_6_nolim_screen) + (1|study_id))
# 
# # df_il_6 = fit_lmer(form, df= mh_il6_vitd, label = "IL-6 (no limits)")
# 
# 
# fit_il6 <- lmer(formula = log(il_6_nolim) ~ rand*time + stratification + 
#                   log(il_6_nolim_screen) + (1|study_id),
#                 data = mh_il6_vitd) 
# 
# 
# 
# 
# # hist(residuals(fit_il6))
# # summary(fit_il6)$coefficients
# # anova(fit_il6)
# (p <- summary(pairs(emmeans(fit_il6, ~ rand * time ), 
#                     simple="rand")))  #contrasts PINK - BLUE i.e. 1000 - 500
# 
# 
# estimate3<-as.numeric(p[1,"estimate"])  #estimate  at 6 weeks
# se3<-as.numeric(p[1,"SE"]) #se at 6 weeks
# p3 <-format.pval(p[1,"p.value"], eps=0.001,2)
# 
# estimate6<-as.numeric(p[2,"estimate"])  #estimate at 12m
# se6<-as.numeric(p[2,"SE"]) #se at 12m
# p6 <-format.pval(p[2,"p.value"], eps=0.001,2)
# 
# time=c("6 weeks", "12 months")
# estimate = rbind(estimate3,estimate6)
# se=rbind(se3,se6)
# p=rbind(p3,p6)
# 
# df_il_6 = data.frame(Outcome = c("IL-6",""),
#                      Contrast = c("500mg vs 1000mg",""), 
#                      Level=time, Estimate = exp(-estimate), 
#                      lower = exp(-estimate - 1.96*se),
#                      upper = exp(-estimate + 1.96*se),
#                      se = se, p=p)
# 
# rownames(df_il_6) <- c()

df_25ohd2 = lmer_nolim(x25ohd2fu, "25(OH)D2", x25ohd2_screen, 
                       mh_il6_vitd, limit=3.0)

df_25ohd3 = lmer_nolim(x25ohd3fu, "25(OH)D3", x25ohd3_screen, 
                       mh_il6_vitd, limit=0)

df_epi25ohd3 = lmer_nolim(epi25ohd3fu, "Epi25(OH)D3", epi25ohd3_screen, 
                          mh_il6_vitd, limit=2.0)
df_il_6 = lmer_nolim(il_6, "IL-6", il_6_screen, mh_il6_vitd, 0.2)


# ---- tobit ----
#https://stats.stackexchange.com/questions/581866/how-to-do-post-hoc-analysis-with-contrasts-of-a-tobit-model-censreg-package
# https://stats.stackexchange.com/questions/624699/post-hoc-analysis-for-tobit-model-using-censreg-in-r

#set up data for tobit regression
#Epi 25(OH)D3 lower limit is 2.0
# 25(OH)D2 lower limit is 3.0
# IL-6 lower limit is 0.2

pmh <- pdata.frame(mh_il6_vitd %>% 
                     arrange(id, time)
                   , c( "id","time" ) )



# summary(mh_il6_vitd$epi25ohd3fu_cens)
# summary(mh_il6_vitd$epi25ohd3_screen_cens)
# hist(log(mh_il6_vitd$epi25ohd3fu_cens))



# 25(OH)D2

form = as.formula(log(x25ohd2fu_cens) ~ rand*time + stratification + log(x25ohd2_screen_cens)) 

df_25ohd2_cens = fit_censreg(form, pmh, leftlim = log(3), 
                             label = "25(OH)D2 (Censored)")



# Epi 25(OH)D3
form = as.formula(log(epi25ohd3fu_cens) ~ rand*time + stratification + log(epi25ohd3_screen_cens) )

df_epi25ohd3_cens = fit_censreg(form, pmh, leftlim = log(2), 
                                label = "Epi25(OH)D3 (Censored)")



# hist(mh_il6_vitd$il_6)

form = as.formula(log(il_6) ~ rand*time + stratification + log(il_6_screen))

df_il_6_cens = fit_censreg(form, pmh, leftlim = log(0.2), 
                           label = "IL-6 (Censored)")



# ---- tab_comb ----
##combine

tab_cens <- bind_rows(df_25ohd2, df_25ohd2_cens, 
                      df_epi25ohd3, df_epi25ohd3_cens,
                      df_25ohd3, 
                      df_il_6,df_il_6_cens) %>%
  mutate(
    est_ci = paste0(round(Estimate, 2), " (",
                    round(lower, 2), ", ",
                    round(upper, 2), ")")) 


# tab_cens %>%
#   select(Outcome, Contrast, Level, est_ci, p) %>%
#   flextable %>%
#   autofit() %>%
#   theme_box() %>%
#   set_header_labels(est_ci = "RR (95% CI)")


