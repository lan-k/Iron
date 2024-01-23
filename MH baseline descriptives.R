##baseline descriptives for MH variables
#BLUE is 500mg
#PINK is 1000 mg

#---- mh_data_base ----

rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(haven)
library(flextable)
library(labelled)

quantile_str <- function(x, probs = c(0.25, 0.5, 0.75)) {
  value = quantile(x, probs, na.rm=T)
  
  med_iqr = paste0(as.character(round(value[2])), " (", 
                   as.character(round(value[1])), "-",
                   as.character(round(value[3])), ")")
         
  return(med_iqr)
}

mh <- read_sas(data_file = "../iron_mh.sas7bdat")

base <- mh %>%
  filter(time==3) %>%
  select(Study_ID, rand, contains("BL"))

names(base)


base_n = base %>%
  group_by(rand) %>%
  summarise(n=n()) %>%
  ungroup()


base_desc = base %>%
  group_by(rand) %>%
  reframe(across(contains("BL_"), quantile_str)) %>%
  ungroup()

base_desc_long <- base_desc %>%
  pivot_longer(cols = !rand,  names_to = "Variable", values_to = "med_iqr") 

base_tab <- full_join(base_desc_long %>% filter(rand=="BLUE") %>%
                        select(!rand) %>%
                        rename(med_iqr_BLUE=med_iqr),
                      base_desc_long %>% filter(rand=="PINK") %>%
                        select(!rand) %>% 
                        rename(med_iqr_PINK=med_iqr),
                      by="Variable")

blue_lab = paste0("500 mg"," (n=", as.character(base_n %>% filter(rand=="BLUE") %>% pull(n)), ")")
pink_lab = paste0("1000 mg"," (n=", as.character(base_n %>% filter(rand=="PINK") %>% pull(n)), ")")

base_tab <- base_tab %>% 
  mutate(Variable=gsub("BL|_|total|MH","", Variable),
         Variable=gsub("DDEP"," - Depression", Variable),
         Variable=gsub("AANX"," - Anxiety", Variable)) 


base_tab %>%
  flextable() %>%
  autofit() %>%
  set_header_labels(Variable = "Variable",
             med_iqr_BLUE = blue_lab,
             med_iqr_PINK = pink_lab) %>%
  set_caption(caption = "Baseline Mental Health Descriptives")

knitr::kable(base_tab)
