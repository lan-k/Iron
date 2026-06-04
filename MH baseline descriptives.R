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


# ---- bio_desc_p ----

# 4/6/2026 add p-values for assessor's comments
rm(list=ls())
library(arsenal)
library(janitor)

contvars = c("Age","BMI", "Education_years", 
             "GA_Screening", "GA_FirstIx",
             "Hb_Screening","Ferritin_Screening","Iron_screening",
             "TF_Screening","TSAT_Screening",
             "MCH_Screening","Serumfolate_Screen2",
             "serumb12_screening","activeb12_screening",
             "CRP_Screening", "il_6_screen",
             "BL_MH_HAMD_DEP","BL_MH_EPDS")
catvars = c("Twins","Anaemia", "income","Supp_code"
)

vars = c(contvars, catvars)

mycontrols  <- tableby.control(test=T, total=F,
                               na.action=na.tableby(F),
                               numeric.test = "wt",
                               cat.test = "chisq",
                               chisq.correct = FALSE,
                               cat.simplify = T,
                               numeric.stats=c("N","Nmiss", "medianq1q3"), #,"min","max"
                               cat.stats=c("Nmiss2","countpct"),
                               stats.labels=list(N="N",Nmiss = "Missing",
                                                 Nmiss2="Missing",medianq1q3='Median (IQR)'),
                               digits=0, digits.p=2, digits.pct=1)



iron <- read.csv(file="../Iron dose RCT_audit_24jan_TIDY.csv", na = "9999") %>%
  select(Study_ID,  Age, BMI, Anaemia, income, 
         Education_years, Education_level, Twins, Supp_code, GA_FirstIx, 
         Serumfolate_Screen2,
         matches("Screening|screening"))



base_desc <- read_sas(data_file = "../iron_mh.sas7bdat") %>%
  filter(time==3) %>%
  select(Study_ID,rand,  contains("BL"))

base_desc = left_join(iron, base_desc) 

load(file="mh_il6_vitd.rds" )



desc = left_join(base_desc, 
                 mh_il6_vitd %>% filter(time == 1) %>%
                   select(study_id, il_6_screen), by=c("Study_ID" = "study_id"))


desc = desc %>%
  rename(#Serumfolate = Serumfolate_screening,
         Ferritin = Ferritin_Screening) %>%
  mutate(#Serumfolate_screening = as.numeric(gsub("\\D", "", Serumfolate)),
         # Serumfolate_screening = ifelse(Serumfolate_screening > 1000, NA, Serumfolate_screening),
         Ferritin_Screening     = as.numeric(gsub("\\D", "", Ferritin))) %>%
  select(!Ferritin)



summary(desc$Serumfolate_Screen2)
summary(desc$Ferritin_Screening)

tabdata = desc %>% 
  select(Study_ID, rand, all_of(vars)) %>%
  mutate(income = case_when(is.na(income) ~ 0, 
                            income == 6 ~ 5, 
                            TRUE ~ income),
         Supp_code = ifelse(is.na(Supp_code), 5,Supp_code),
         
         across(all_of(catvars), ~factor(.)))


formula <- as.formula(paste0("rand ~  ", paste0(vars, collapse="+")))

tabdesc <- tableby(formula, data=tabdata , 
                       control=mycontrols)


summary(tabdesc, text = T, pfootnote = T)

summary(padjust(tabdesc, method= "bonf",
                suffix = " (adjusted for multiple comparisons)"), 
        title='Patient Characteristics by group', 
        text = T, pfootnote = T)



tabdesc2 <- tableby(rand ~ Hb_Screening + log(il_6_screen), data=tabdata, test=T)
summary(tabdesc2, text = T, pfootnote = T)

                      