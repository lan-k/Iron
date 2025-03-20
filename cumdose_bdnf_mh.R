###biomarker data
##topup_orig: 0=no; 1= yes; 2= required but not administered; 3 = administered by not required
#BLUE is 500mg
#PINK is 1000 mg
#prev_topup is whether a topup was given at the pervious visit
#ntopup is number of topup at current visit
# prev_ntopup is number of topup at previous visit
#cumdose is cumulative dose at previous visit

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


 
# rct <- read.csv("../../biomarkers/Iron dose RCT_audit_24jan_TIDY.csv",
#                 stringsAsFactors = F, na.strings="") %>% 
#   janitor::clean_names() 

iron <- read_sas(data_file = "../iron_long.sas7bdat") %>%
  mutate(topup_given = ifelse(topup_orig %in% c(1,3), 1, 0),
         stratification = factor(stratification)) %>%
  group_by(Study_ID) %>%
  mutate(ntopup = cumsum(topup_given),
         prev_ntopup =  lag(ntopup),
         prev_ntopup = ifelse(is.na(prev_ntopup), 0, prev_ntopup),
         iron_cum = ifelse(rand=="BLUE", ntopup + 1, ntopup  + 2 ),
         cumdose = ifelse(row_number() == 1, iron_cum*500,lag(iron_cum)*500)) %>%
  ungroup()
#cumdose is cumulative dose BEFORE time point


# dose <- read.csv(file="../../biomarkers/PO_cumulative_dose_days.csv")
# 
# dose_long <- dose %>%
#   pivot_longer(cols=contains("time"),names_to = "time", names_transform = readr::parse_number,
#                values_to="time_since_last_dose") %>%
#   select(!contains("cumdose|PO"),  !Anemia)

##these 
# dose_long <- bind_rows(dose %>% select(Study_ID, cumdose_FU1, time_last_dose_FU1) %>%
#                          mutate(time = 3) %>%
#                          rename(cumdose=cumdose_FU1, time_last_dose = time_last_dose_FU1),
#                        dose %>% select(Study_ID, cumdose_FU2, time_last_dose_FU2) %>%
#                          mutate(time = 6) %>%
#                          rename(cumdose=cumdose_FU2, time_last_dose = time_last_dose_FU2)) %>%
#   arrange(Study_ID, time)

# dose_dates <- read.csv(file="../../biomarkers/variable_dates_13_02_2023.csv", na.string="") %>%
#   select(Study_ID, contains("_date"), contains("_FU")) %>%
#   filter(!is.na(Study_ID))
# 
# #convert to long
# 
# dose_dates_long <- dose_dates %>%
#   pivot_longer(!Study_ID, values_to = "Date", names_to = "time")


  
last_dose <- iron %>%
  filter(topup_given == 1 | time %in% c(3,6)) %>%
  select(Study_ID,time, rand,stratification, topup_given,iron_cum,  cumdose) %>%
  group_by(Study_ID) %>%
  mutate(prev_topup= lag(topup_given),
         prev_meas = lag(time),
         prev_topup_time = ifelse(prev_topup == 1, prev_meas, lag(prev_meas)))


iron_topup <- left_join(iron, last_dose %>% 
                          select(!c(prev_meas))) %>%
  # left_join(dose_long) %>%
  filter(time %in% c(3,6)) %>%
  mutate(prev_topup_time = factor(ifelse(is.na(prev_topup_time), 0, prev_topup_time), 
                                  levels=0:5, ordered = TRUE))
         
table(iron_topup$time, iron_topup$prev_topup_time)

#if no topup given then prev time is at randomisation


mh <- read_sas(data_file = "../iron_mh.sas7bdat") %>%
  filter(time %in% c(3,6)) %>%
  mutate(time=ifelse(time==3, "FU1","FU2"))

ids <- mh %>% pull(Study_ID) %>% unique()

agp <- read.csv("../../biomarkers/04_agp.csv",stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names()
cortisol <- read.csv("../../biomarkers/10_cortisol.csv",stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names()


bdnf <- read.csv("../../biomarkers/19_BDNF_complete_dataset_replaced_15_BDNF_in_R.csv",
                 stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names() %>%
  mutate(bdnf_pg_ml = bdnf_pg_ml/1000)

dat_bdnf <- read.csv("../../biomarkers/joined_left_1.csv",stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names() %>%
  mutate(randomisation = relevel(factor(randomisation), ref="PINK")) %>%
  filter(study_no %in% ids)


dat_bdnf <- dat_bdnf %>%
  left_join(bdnf) %>%
  left_join(cortisol) %>%
  rename(time = sample) %>%
  select(!c(contains("study_no2"), contains("sample"))) %>%
  filter(!is.na(randomisation)) %>%
  mutate(rand=relevel(factor(randomisation), ref="PINK"))
  


base = dat_bdnf %>% filter(time == "SCR")

dat <- dat_bdnf %>% filter( time != "SCR") %>%
  left_join(base %>% select(study_no, bdnf_pg_ml) %>%
              rename(scr_bdnf = bdnf_pg_ml), by = "study_no") %>%
  left_join(iron_topup %>% select(Study_ID,time,stratification, 
                                  # time_since_last_dose,
                                  ntopup,prev_ntopup,
                                  topup_given, cumdose, prev_topup,
                                  prev_topup_time) %>%

              mutate(time=ifelse(time==3, "FU1","FU2"),
                     prev_topup = ifelse(is.na(prev_topup), 0, prev_topup),
                     prev_topup = relevel(factor(prev_topup), ref="0")), 
            by=c("study_no"="Study_ID", "time")) %>%
  left_join(mh %>% select(Study_ID,time, Age, income, Education_years), 
            by=c("study_no"="Study_ID", "time") ) 



cord <- dat %>%
  filter(time == "CB") %>%
  select(!stratification) %>%
  left_join(iron %>% select(Study_ID, stratification) %>% unique(), 
            by=c("study_no"="Study_ID"))

dat_noCB <- dat %>% filter(time != "CB") %>%
  mutate(time=relevel(factor(time), ref="FU1"),
         bdnf_change = bdnf_pg_ml - scr_bdnf)


table(dat_noCB$time)

tapply(dat_noCB$bdnf_pg_ml, dat_noCB$randomisation, summary)

tapply(dat_noCB$bdnf_change, dat_noCB$randomisation, summary)
tapply(dat_noCB$bdnf_change, dat_noCB$time, summary)


mh <- left_join(mh, dat_noCB %>% select(study_no,time, cumdose, ntopup ),
                by=c("Study_ID"="study_no","time")) 
#---- desc_bdnf ----

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



bdnf_desc = dat_bdnf %>%
  group_by(rand, time) %>%
  reframe(across(contains("bdnf"), \(x) quantile_str(x, digits=1))) %>%
  ungroup()


bdnf_desc_wide <- bdnf_desc %>%
  pivot_wider(names_from = "rand", values_from = "bdnf_pg_ml") 


bdnf_desc_wide %>%
  flextable() %>%
  autofit() %>%
  set_header_labels(time = "Time point",
                    BLUE = blue_lab,
                    PINK = pink_lab) %>%
  set_caption(caption = "BDNF (pg/ml) (N)")

knitr::kable(bdnf_desc_wide)

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
  
  form = formula(paste0(outvar, " ~",  exposure,"  + time + " , paste(adjvars, collapse = "+"),
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


