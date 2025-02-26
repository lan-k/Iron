###biomarker data

#BLUE is 500mg
#PINK is 1000 mg
#mediation analysis of IV dose on MH_HAMD_DEP and MH_EPDS via Vitamin D and IL-6 levels
# note: MH_EPDS_12 is from 12 week visit (not 12 month)
# exposure variable is treatment group

#limits of detection: Epi25(OH)d3 <2.0, 25(OH)d2 <3.0, 25(OH)d3 no limit in data
#IL-6 there seem to be 2 limits: 0.9 and 0.2

#Model 1: biomarker as the outcome and IVI as the exposure
#Model 2: mediation analysis. MH as the outcome, biomarkers as the mediator, ICI as the exposure

#Model 1: Tobit regression for left truncated biomarkers 
# + lmer for normally distributed + estimation below the threshold using random number, compare

#Model 2: use mediate R package for 25(OH)d3
# see results of model 1 for other biomarkers, use random estimation under threshold if required


#---- biomarker_data ----

rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(haven)
library(flextable)



mh <- read_sas(data_file = "../iron_mh.sas7bdat") %>%
  filter(time %in% c(3,6)) %>%
  mutate(time=ifelse(time==3, "FU1","FU2")) %>% 
  janitor::clean_names()


hist(mh$mh_epds) #change from baseline
hist(mh$mh_hamd_dep)  #change from baseline

ids <- mh %>% pull(study_id) %>% unique()



il6 <- read.csv("../../VitD and IL6/18_inflammatory_panel.csv",
                stringsAsFactors = F, na.strings="") %>% 
  janitor::clean_names() %>%
  mutate(il_6 = case_when(il_6_pg_ml == "<=0.9" ~ 0.5 + 0.4*runif(1),
                          il_6_pg_ml == "<=0.2" ~ 0.2,
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

mh_il6_vitd <- mh_il6_vitd %>%
  mutate(FU = time,
         time = factor(ifelse(time=="FU1", 1, 2)))

save(base,mh_il6_vitd, file="mh_il6_vitd.rds" )

