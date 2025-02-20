###biomarker data
##topup_orig: 0=no; 1= yes; 2= required but not administered; 3 = administered by not required
#BLUE is 500mg
#PINK is 1000 mg
#mediation analysis of IV dose on MH_HAMD_DEP and MH_EPDS via Vitamin D and IL-6 levels
# note: MH_EPDS_12 is from 12 week visit (not 12 month)
# exposure variable is treatment group


#---- biomarker_data ----

rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(haven)
library(lme4)
library(lmerTest)
library(emmeans)
library(splines)
library(flextable)



mh <- read_sas(data_file = "../iron_mh.sas7bdat") %>%
  filter(time %in% c(3,6)) %>%
  mutate(time=ifelse(time==3, "FU1","FU2")) %>% 
  janitor::clean_names()

ids <- mh %>% pull(study_id) %>% unique()



il6 <- read.csv("../../VitD and IL6/18_inflammatory_panel.csv",
                stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names() %>%
  mutate(il_6 = case_when(il_6_pg_ml == "<=0.9" ~ 0.89,
                          il_6_pg_ml == "<=0.2" ~ 0.19,
                          TRUE ~ as.numeric(il_6_pg_ml))) %>%
  rename(time = sample, study_id = study_no) %>% 
  select(study_id, time, il_6) %>%
  arrange(study_id, time) 

il6_base <- il6 %>% filter(time == "SCR") %>%
  unique() %>%
  rename(il_6_screen = il_6) %>%
  select(!time)

il6_cb <- il6 %>% filter(time == "CB") %>%
  unique()

il6_fu <- il6 %>% filter(time %in% c("FU1", "FU2"))   %>%
  unique()  


vitd <- read.csv("../../VitD and IL6/PO_FGF_Ca_Haem_vitD_long_for_stats_step2.csv",
                 stringsAsFactors = F, na.strings="NA") %>% 
  janitor::clean_names() %>%
  arrange(study_id, fu) %>%
  mutate(time=ifelse(fu == 1, "FU1", "FU2")) 



mh_il6_vitd <- mh %>%
  select(study_id, rand, time, stratification, age, income, education_years, 
         matches("epds|hamd") ) %>%
  left_join(il6_fu ) %>%
  left_join(il6_base) %>%
  left_join(vitd %>% select(study_id, time, contains("oh"))) %>%
  arrange(study_id, time) %>%
  filter(!is.na(rand)) %>%
  mutate(rand=relevel(factor(rand), ref="PINK"))



base = mh_il6_vitd %>%
  select(study_id, rand, contains("screen")) %>%
  unique()

colnames(base) <-  sub("_screen.*", "", colnames(base))

#---- desc ----



quantile_str <- function(x, probs = c(0.25, 0.5, 0.75), digits=0) {
  value = quantile(x, probs, na.rm=T, digits=digits)
  n=sum(!is.na(x))
  
  med_iqr = paste0(as.character(round(value[2],digits=digits)), " (", 
                   as.character(round(value[1],digits=digits)), "-",
                   as.character(round(value[3],digits=digits)), ")",
                   " (N=",n,")")
  
  return(med_iqr)
}

base_n = base %>%
  group_by(rand) %>%
  summarise(n=n()) %>%
  ungroup()


blue_lab = paste0("500 mg"," (n=", as.character(base_n %>% filter(rand=="BLUE") %>% pull(n)), ")")
pink_lab = paste0("1000 mg"," (n=", as.character(base_n %>% filter(rand=="PINK") %>% pull(n)), ")")


base_desc <- base %>%
  group_by(rand) %>%
  reframe(across(matches("il_6|oh"), \(x) quantile_str(x, digits=2))) %>%
  ungroup() %>%
  mutate(time="SCR")


base_desc_wide <- base_desc %>%
  pivot_wider(names_from = "rand", 
              values_from = c("il_6","epi25ohd3","x25ohd2", "x25ohd3") )


desc = mh_il6_vitd %>%
  select(!contains("screen")) %>%
  group_by(rand, time) %>%
  reframe(across(matches("il_6|oh"), \(x) quantile_str(x, digits=2))) %>%
  ungroup()

colnames(desc) <-  sub("fu", "", colnames(desc))


desc_wide <- desc %>%
  pivot_wider(names_from = "rand", 
              values_from = c("il_6","epi25ohd3","x25ohd2", "x25ohd3") )

desc_wide <-bind_rows(base_desc_wide, desc_wide) 



desc_wide %>%
  flextable() %>%
  autofit() %>%
  set_header_labels(time = "Time point",
                    BLUE = blue_lab,
                    PINK = pink_lab) %>%
  set_caption(caption = "Biomarkers (N)")

knitr::kable(desc_wide)

#---- desc_cumdose ----

desc_cumdose <- dat_noCB %>%
  group_by(rand, time) %>%
  summarise(cumdose = mean(cumdose, na.rm=T),
            ntopup = mean(ntopup, na.rm=T),
            prev_ntopup = mean(prev_ntopup, na.rm=T)) %>%
  ungroup() 
##sd not working

pink <- dat_noCB %>% filter(rand == "PINK")
blue <- dat_noCB %>% filter(rand == "BLUE")

tapply(pink$cumdose, pink$time, sd)
tapply(blue$cumdose, blue$time, sd)


tapply(pink$ntopup, pink$time, sd)
tapply(blue$ntopup, blue$time, sd)

tapply(pink$prev_ntopup, pink$time, sd)
tapply(blue$prev_ntopup, blue$time, sd)

#---- lmer_bdnf ----

dat_noCB <- dat_noCB %>% 
  mutate(cumdosesm = cumdose/500) 

#randomisation group
##cord blood
fitrand_cb <- lm(log(bdnf_pg_ml) ~ rand + stratification + log(scr_bdnf),
                 data = cord)

(s_cb <- summary(fitrand_cb)$coefficients)

df_cb <- data.frame(Covariate = "Cord blood: 500mg vs 1000mg",
                      Estimate = s_cb["randBLUE","Estimate"],
                      se = s_cb["randBLUE","Std. Error"],
                      p=format.pval(s_cb["randBLUE","Pr(>|t|)"], eps=0.001,2))



#maternal
fitrand <- lmer(formula = bdnf_change ~ rand*time + stratification + 
                  log(scr_bdnf) + (1|study_no),
                data = dat_noCB) 
summary(fitrand)$coefficients
anova(fitrand)
(p <- summary(pairs(emmeans(fitrand, ~ rand * time ), 
                    simple="rand")))  #contrasts

estimate3<-as.numeric(p[1,"estimate"])  #estimate  at 6 weeks
se3<-as.numeric(p[1,"SE"]) #se at 6 weeks
p3 <-format.pval(p[1,"p.value"], eps=0.001,2)

estimate6<-as.numeric(p[2,"estimate"])  #estimate at 12m
se6<-as.numeric(p[2,"SE"]) #se at 12m
p6 <-format.pval(p[2,"p.value"], eps=0.001,2)

time=c("6 weeks", "12 months")
estimate = rbind(estimate3,estimate6)
se=rbind(se3,se6)
p=rbind(p3,p6)

df_rand = data.frame(Covariate = c("Maternal: 500mg-1000mg",""), 
                     Level=time, Estimate = -1*estimate, 
                se = se, p=p)



##cumulative dose

fitdose <- lmer(formula = bdnf_change ~ cumdosesm + stratification + 
                  log(scr_bdnf) + (1|study_no),
                data = dat_noCB) 

(s_dose <- summary(fitdose)$coefficients)
anova(fitdose)
hist(residuals(fitdose))

df_dose <- data.frame(Covariate = c("Cumulative dose (increase in 500mg)"),
                      Estimate = s_dose["cumdosesm","Estimate"],
                      se = s_dose["cumdosesm","Std. Error"],
                      p=format.pval(s_dose["cumdosesm","Pr(>|t|)"], eps=0.001,2))


##number of topups
fittopup <- lmer(formula = bdnf_change ~ prev_ntopup + stratification + 
                  log(scr_bdnf) + (1|study_no),
                data = dat_noCB) 

(s_topup <- summary(fittopup)$coefficients)
anova(fittopup)
hist(residuals(fittopup))

df_topup <- data.frame(Covariate = "Number of repeat infusions",
                    Estimate = s_topup["prev_ntopup","Estimate"],
                    se = s_topup["prev_ntopup","Std. Error"],
                    p=format.pval(s_topup["prev_ntopup","Pr(>|t|)"], eps=0.001,2))


##randomisation group and previous topup

fitrand_prev <- lmer(formula = bdnf_change ~ rand*prev_topup + stratification + 
                  log(scr_bdnf) + (1|study_no),
                data = dat_noCB) 

summary(fitrand_prev)


anova(fitrand_prev)
(pprev <- summary(pairs(emmeans(fitrand_prev, ~ rand * prev_topup ), 
                    simple="prev_topup")))  #contrasts


estimate1000<-as.numeric(pprev[1,"estimate"])  
se1000<-as.numeric(pprev[1,"SE"]) 
p1000 <-format.pval(pprev[1,"p.value"], eps=0.001,2)

estimate500<-as.numeric(pprev[2,"estimate"])  
se500<-as.numeric(pprev[2,"SE"]) 
p500 <-format.pval(pprev[2,"p.value"], eps=0.001,2)

level=c("500 mg", "1000 mg")
estimate = rbind(estimate500,estimate1000)
se=rbind(se500,se1000)
p=rbind(p500,p1000)



df_prev = data.frame(Covariate = c("Previous top-up",""), Level=level, Estimate = -1*estimate, 
                     se = se, p=p)


##previous dose and time since last dose

# fitrand_prev_time <- lmer(formula = bdnf_change ~ randomisation*prev_topup + stratification + 
#                       bs(time_last_dose) + log(scr_bdnf) + (1|study_no),
#                      data = dat_noCB) 
# 
# summary(fitrand_prev_time)
# anova(fitrand_prev_time)
# (pprev <- summary(pairs(emmeans(fitrand_prev_time, ~ randomisation * prev_topup ), 
#                         simple="prev_topup")))


##combine

tab_bdnf <- bind_rows(df_cb,df_rand, df_dose, df_topup, df_prev) %>%
  mutate(lower = Estimate - 1.96*se,
         upper = Estimate + 1.96*se,
         est_ci = paste0(round(Estimate, 2), " (",
                         round(lower, 2), ", ",
                         round(upper, 2), ")"),
         est_ci=ifelse(grepl("Cord",Covariate),
                       paste0(round(exp(Estimate), 2), " (",
                              round(exp(lower), 2), ", ",
                              round(exp(upper), 2), ")"),
                       est_ci)) 


tab_bdnf %>%
  select(Covariate, Level, est_ci, p) %>%
  flextable %>%
  autofit() %>%
  theme_box() %>%
  set_header_labels(est_ci = "Estimate (95% CI)")


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


