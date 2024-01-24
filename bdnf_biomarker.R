###biomarker data
##topup_orig: 0=no; 1= yes; 2= required but not administered; 3 = administered by not required
#BLUE is 500mg
#PINK is 1000 mg

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

 
# rct <- read.csv("../../biomarkers/Iron dose RCT_audit_24jan_TIDY.csv",
#                 stringsAsFactors = F, na.strings="") %>% 
#   janitor::clean_names() 

iron <- read_sas(data_file = "../iron_long.sas7bdat") %>%
  mutate(topup_given = ifelse(topup_orig %in% c(1,3), 1, 0),
         stratification = factor(stratification)) %>%
  group_by(Study_ID) %>%
  mutate(iron_cum = cumsum(topup_given),
         iron_cum = ifelse(rand=="BLUE", iron_cum + 1, iron_cum  + 2 ),
         cumdose = 500*iron_cum) %>%
  ungroup()

time_since_last_dose <- iron %>%
  filter(topup_given == 1 | time %in% c(3,6)) %>%
  select(Study_ID,time, rand,stratification, topup_given,iron_cum,  cumdose) %>%
  group_by(Study_ID) %>%
  mutate(prev_topup= lag(topup_given),
         prev_meas = lag(time),
         prev_topup_time = ifelse(prev_topup == 1, prev_meas, lag(prev_meas)))


iron_topup <- left_join(iron, time_since_last_dose %>% 
                          select(!c(prev_meas))) %>%
  filter(time %in% c(3,6)) %>%
  mutate(prev_topup_time = factor(ifelse(is.na(prev_topup_time), 0, prev_topup_time), 
                                  levels=0:5, ordered = TRUE))
         
table(iron_topup$time, iron_topup$prev_topup_time)

#if no topup given then prev time is at randomisation


mh <- read_sas(data_file = "../iron_mh.sas7bdat") %>%
  filter(time %in% c(3,6)) %>%
  mutate(time=ifelse(time==3, "FU1","FU2"))

agp <- read.csv("../../biomarkers/04_agp.csv",stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names()
cortisol <- read.csv("../../biomarkers/10_cortisol.csv",stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names()
bdnf <- read.csv("../../biomarkers/19_BDNF_complete_dataset_replaced_15_BDNF_in_R.csv",
                 stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names() %>%
  mutate(bdnf_pg_ml = bdnf_pg_ml/1000)

dat <- read.csv("../../biomarkers/joined_left_1.csv",stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names() %>%
  mutate(randomisation = relevel(factor(randomisation), ref="BLUE"))


dat <- dat %>%
  left_join(bdnf) %>%
  left_join(cortisol) %>%
  rename(time = sample) %>%
  select(!c(contains("study_no2"), contains("sample"))) %>%
  filter(!is.na(randomisation))
  


base = dat %>% filter(time == "SCR")

dat <- dat %>% filter( time != "SCR") %>%
  left_join(base %>% select(study_no, bdnf_pg_ml) %>%
              rename(scr_bdnf = bdnf_pg_ml), by = "study_no") %>%
  left_join(iron_topup %>% select(Study_ID,time,stratification, 
                                  topup_given, cumdose, prev_topup,
                                  prev_topup_time) %>%

              mutate(time=ifelse(time==3, "FU1","FU2"),
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


#---- lmer_bdnf ----
##cord blood
fitrand_cb <- lm(log(bdnf_pg_ml) ~ randomisation + stratification + log(scr_bdnf),
                 data = cord)

summary(fitrand_cb)


fitrand <- lmer(formula = bdnf_change ~ randomisation*time + stratification + 
                  log(scr_bdnf) + (1|study_no),
             data = dat_noCB) 
summary(fitrand)
anova(fitrand)
(p <- summary(pairs(emmeans(fitrand, ~ randomisation * time ), 
                            simple="randomisation")))  #contrasts


##cumulative dose
fitdose <- lmer(formula = bdnf_change ~ log(cumdose) + stratification + 
                  log(scr_bdnf) + (1|study_no),
                data = dat_noCB) 

summary(fitdose)
anova(fitdose)


##randomisation group and previous topup

fitrand_prev <- lmer(formula = bdnf_change ~ randomisation*prev_topup + stratification + 
                  log(scr_bdnf) + (1|study_no),
                data = dat_noCB) 

summary(fitrand_prev)
anova(fitrand_prev)
(pprev <- summary(pairs(emmeans(fitrand_prev, ~ randomisation * prev_topup ), 
                    simple="prev_topup")))  #contrasts





