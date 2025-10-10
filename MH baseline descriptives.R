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

quantile_str_n <- function(x, probs = c(0.25, 0.5, 0.75), digits=2) {
  value = quantile(x, probs, na.rm=T)
  n=sum(!is.na(x))
  
  med_iqr = paste0(as.character(round(value[2], digits=digits)), " (", 
                   as.character(round(value[1], digits=digits)), "-",
                   as.character(round(value[3], digits=digits)), ") (n=",
                   n, ")")
         
  return(med_iqr)
}


quantile_str <- function(x, probs = c(0.25, 0.5, 0.75)) {
  value = quantile(x, probs, na.rm=T)
  
  med_iqr = paste0(as.character(round(value[2])), " (", 
                   as.character(round(value[1])), "-",
                   as.character(round(value[3])), ") ")
  
  return(med_iqr)
}


# ---- mh_desc ----
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
  reframe(across(contains("BL_"), quantile_str_n)) %>%
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

# ---- bio_desc ----

##median IQr of Vit D and Il-D

load(file="mh_il6_vitd.rds" )

bio_base_desc <- base %>%
  group_by(rand) %>%
  reframe(across(matches("il_6|oh"), quantile_str_n)) %>% #median_iqr
  ungroup() 

bio_base_desc_long <- bio_base_desc %>%
  pivot_longer(cols = !rand,  names_to = "Variable", values_to = "med_iqr") 

bio_base_tab <- full_join(bio_base_desc_long %>% filter(rand=="BLUE") %>%
                        select(!rand) %>%
                        rename(med_iqr_BLUE=med_iqr),
                        bio_base_desc_long %>% filter(rand=="PINK") %>%
                        select(!rand) %>% 
                        rename(med_iqr_PINK=med_iqr),
                      by="Variable")

blue_lab = paste0("500 mg"," (n=", as.character(base_n %>% filter(rand=="BLUE") %>% pull(n)), ")")
pink_lab = paste0("1000 mg"," (n=", as.character(base_n %>% filter(rand=="PINK") %>% pull(n)), ")")

bio_base_tab <- bio_base_tab %>% 
  mutate(Variable=str_to_upper(gsub("x","", Variable)))


bio_base_tab %>%
  flextable() %>%
  autofit() %>%
  set_header_labels(Variable = "Variable",
                    med_iqr_BLUE = blue_lab,
                    med_iqr_PINK = pink_lab) %>%
  set_caption(caption = "Baseline Biomarker Descriptives")

knitr::kable(base_tab)



bio_base_desc_wide <- bio_base_desc %>%
  pivot_wider(names_from = "rand", 
              values_from = c("il_6","epi25ohd3","x25ohd2", "x25ohd3") )




bio_desc_wide %>%
  select(FU, il_6_BLUE,il_6_PINK,epi25ohd3_BLUE, epi25ohd3_PINK,
         x25ohd2_BLUE, x25ohd2_PINK, x25ohd3_BLUE, x25ohd3_PINK) %>%
  flextable() %>%
  autofit() %>%
  set_header_labels(FU = "Time point",
                    il_6_BLUE = "IL-6, \n500 mg",
                    il_6_PINK = "IL-6, \n1000 mg",
                    epi25ohd3_BLUE = "Epi-25(OH)D3, \n500mg",
                    epi25ohd3_PINK = "Epi-25(OH)D3, \n1000mg",
                    x25ohd2_BLUE = "25(OH)D2, \n500mg",
                    x25ohd2_PINK = "25(OH)D2, \n1000mg",
                    x25ohd3_BLUE = "25(OH)D3, \n500mg",
                    x25ohd3_PINK = "25(OH)D3, \n1000mg"
                    ) %>%
  set_caption(caption = "Biomarkers (mean (SD))")

