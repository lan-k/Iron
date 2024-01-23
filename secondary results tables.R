##tables for secondary outcomes
rm(list=ls())
library(dplyr)
library(stringr)
library(tidyr)
library(haven)
library(naniar)
library(scales)



iron <- read_sas("../Data/iron.sas7bdat") %>%
  mutate(rand=ifelse(rand == "BLUE", "500mg", "1000mg"),
         Delivery = case_when(Mode_Delivery == 1 ~ "Normal vaginal",
                              Mode_Delivery == 2 ~ "Instrumental vaginal", 
                              Mode_Delivery == 3 ~ "Episiotomy/vaginal", 
                              Mode_Delivery == 4 ~ "Emergency Caesarean section", 
                              Mode_Delivery == 5 ~ "Elective Caesarean section"),
         PPROM = ifelse(PPROM == 9999, NA, PPROM)) %>%
  replace_with_na(replace = list(9999))

iron_long <- read_sas("../Data/iron_long.sas7bdat") %>%
  mutate(rand=ifelse(rand == "BLUE", "500mg", "1000mg")) %>%
  replace_with_na(replace = list(9999))

iron_blood <- read_sas("../Data/iron_blood.sas7bdat") %>%
  mutate(rand=ifelse(rand == "BLUE", "500mg", "1000mg")) %>%
  replace_with_na(replace = list(9999))


ferritin_long <- read_sas("../Data/ferritin_long.sas7bdat") %>%
  mutate(rand=ifelse(rand == "BLUE", "500mg", "1000mg")) %>%
  replace_with_na(replace = list(9999))


infant <- read_sas("../Data/infant.sas7bdat") %>%
  mutate(rand=ifelse(rand == "BLUE", "500mg", "1000mg")) %>%
  replace_with_na(replace = list(9999))



# function for rounding numbers with zeros kept
roundz = function(x, digits){
  dformat = paste('%.', digits, 'f', sep='')
  x = sprintf(dformat, round(x, digits))
  return(x)
}


##fisher test for PPROM
table(iron$rand, iron$PPROM)
fisher.test(factor(iron$rand), factor(iron$PPROM)) #p=0.004

#### secondary results ####
bloodtab <- function(outvar, df, var) {
  
  out <- iron_sec %>% 
    filter(Outcome == outvar)
  
  obs <- df %>%
    select(time, rand, {{var}}) %>%
    group_by(time ,rand) %>%
    filter(!is.na({{var}})) %>%
    summarise(N=n(),mean=mean({{var}}, na.rm=T), 
              sd=sd({{var}}, na.rm=T)) %>%
    mutate(obs_sd = paste0(roundz(mean,digits = 2), " (", 
                           roundz(sd,digits = 2),")",
                           " (",N,")"),
           time=case_when(time == 1 ~ "4 weeks post-randomisation",
                          time == 3 ~ "6 weeks post-partum",
                          time == 4 ~ "3 months post-partum",
                          time == 5 ~ "6 months post-partum",
                          time == 6 ~ "12 months post-partum"))
  
  obs500 <- obs %>%
    filter(rand=="500mg") %>%
    rename(obs_sd_500=obs_sd) %>%
    select(time, obs_sd_500)
  
  obs1000 <- obs %>%
    filter(rand=="1000mg") %>%
    rename(obs_sd_1000=obs_sd) %>%
    select(time, obs_sd_1000)
  
  obs_sec <- out %>% select(Outcome, time) %>%
    left_join(obs500) %>%
    left_join(obs1000) %>%
    left_join(out %>% 
                select(Estimate,p, time)) %>%
    mutate(Outcome = ifelse(time == "4 weeks post-randomisation", 
                            Outcome, NA_character_))
  
  return(obs_sec)
  
}




# iron_obs <- read.csv(file="U:/NALHN/Iron study/Results/secondary/observed iron blood levels 15FEB22.csv",  
#                      stringsAsFactors = F) %>%
#   select(!X_TYPE_) %>%
#   rename(N=X_FREQ_)

iron_sec <- read.csv(file="../Results/secondary/iron blood group diff mixed 15FEB22.csv",  
                     stringsAsFactors = F) 


iron_sec <- iron_sec %>%
  filter(grepl("[0-9]", Label)) %>%
  mutate(Outcome = str_to_sentence(word(Label, 1)),
         Outcome = case_when(grepl("0.5",Label) ~ "Ferritin (half threshold)",
                             grepl("rand",Label) ~ "Ferritin (random method)",
                             Outcome == "Tsat" ~ "TSAT",
                             Outcome == "Tf" ~ "Transferrin",
                             TRUE ~ Outcome),
         Est =paste0(roundz(Estimate,digits = 2), " (", 
                     roundz(Lower,digits = 2),", ",
                     roundz(Upper,digits = 2),")"),
         time=case_when(grepl("4wee",Label) ~ "4 weeks post-randomisation",
                        grepl("6 w",Label) ~ "6 weeks post-partum",
                        grepl("3 mo",Label) ~ "3 months post-partum",
                        grepl("6 mo",Label) ~ "6 months post-partum",
                        grepl("12 m",Label) ~ "12 months post-partum")) %>%
  select(Outcome, time, Est, Probt) %>%
  filter(!is.na(time)) %>%
  rename(Estimate=Est, p=Probt)


irontab <- bloodtab("Iron", iron_blood, iron)
tftab <- bloodtab("Transferrin", iron_blood, tf)
tsattab <- bloodtab("TSAT", iron_blood, tsat)
ferrhalftab <- bloodtab("Ferritin (half threshold)", ferritin_long, ferritin)
ferrrandtab <- bloodtab("Ferritin (random method)", ferritin_long, ferritin) 

iron_sectab <- bind_rows(irontab, tftab, tsattab, ferrhalftab, ferrrandtab)



obstab <- function(outdf, outvar, obsdf, var, cat=1) {
  ##birth and infant outcomes with one time point  
  out <- outdf %>% 
    filter(outcome == outvar)
  
  if (cat != 1) {
    obs <- obsdf %>%
      select(rand, {{var}}) %>%
      mutate(var={{var}},
             var=ifelse(var == 9999, NA_real_,var )) %>%
      group_by(rand) %>%
      filter(var != 9999 | !is.na(var)) %>%
      summarise(N=n(),mean=mean(var, na.rm=T), 
                sd=sd(var, na.rm=T)) %>%
      mutate(obs_sd = paste0(roundz(mean,digits = 1), " (", 
                             roundz(sd,digits = 1),")",
                             " (",N,")"))
    
  } else {
    obs <- obsdf %>%
      select(rand, {{var}}) %>%
      mutate(var={{var}},
             var=ifelse(var == 9999, NA_real_,var )) %>%
      group_by(rand) %>%
      filter(var != 9999 | !is.na(var)) %>%
      summarise(N=n(),sum=sum(var)) %>%
      mutate(prop=sum/N,
             obs_sd = paste0(sum, "/",N, " (", 
                             percent(prop, accuracy = 1),")"))
  }
  
  
  
  obs500 <- obs %>%
    filter(rand=="500mg") %>%
    rename(obs_sd_500=obs_sd) %>%
    select(obs_sd_500)
  
  obs1000 <- obs %>%
    filter(rand=="1000mg") %>%
    rename(obs_sd_1000=obs_sd) %>%
    select(obs_sd_1000)
  
  obs_sec <- out %>% select(outcome) 
  obs_sec$obs_sd_500 = obs500$obs_sd_500
  obs_sec$obs_sd_1000 = obs1000$obs_sd_1000
  obs_sec$Estimate = out$Estimate
  obs_sec$p=out$p
  
  return(obs_sec)
  
}


ordtab <- function(outdf, outvar, obsdf, var, cat=1) {
  ##birth and infant outcomes with one time point  
  out <- outdf %>% 
    filter(outcome == outvar)
  
  
  obs <- obsdf %>%
    select(rand, {{var}}) %>%
    mutate(var={{var}},
           var=ifelse(var == 9999, NA_real_,var )) %>%
    filter(var != 9999 | !is.na(var)) %>%
    group_by(rand, var) %>%
    summarise(n=n()) %>%
    mutate(N=sum(n),
           prop=n/N,
           obs_sd = paste0(n, "/",N, " (", 
                           percent(prop, accuracy = 1),")")) %>%
  ungroup()

  
  obs500 <- obs %>%
    filter(rand=="500mg") %>%
    rename(obs_sd_500=obs_sd) %>%
    select(var, obs_sd_500) %>%
    rename(Response = var)
  
  obs1000 <- obs %>%
    filter(rand=="1000mg") %>%
    rename(obs_sd_1000=obs_sd) %>%
    select(var, obs_sd_1000) %>%
    rename(Response = var)
  
  obs_sec <- out %>% select(outcome, Response) %>%
    full_join(obs500, by="Response")  %>%
    left_join(obs1000, by="Response")  %>%
    left_join(out %>% select(Response, Estimate, p), by="Response")

  
  return(obs_sec)
  
}



### birth outcomes
birth_cont <- read.csv(file="../Results/secondary/Birth Continuous 11FEB22.csv",  
                       stringsAsFactors = F) %>%
  mutate(Estimate=paste0(roundz(as.numeric(MeanEstimate),digits = 3), " (", 
                         roundz(as.numeric(MeanLowerCL),digits = 3),", ",
                         roundz(as.numeric(MeanUpperCL),digits = 3),")")) %>%
  rename(p=ProbChiSq) %>%
  mutate(outcome=gsub("BLUE","500mg",outcome),
         outcome=gsub("PINK","1000mg",outcome))


birth_or <- read.csv(file="../Results/secondary/Birth Odds Ratios 11FEB22.csv",  
                       stringsAsFactors = F) %>%
  mutate(Estimate=paste0(roundz(as.numeric(OddsRatioEst),digits = 3), " (", 
                         roundz(as.numeric(LowerCL),digits = 3),", ",
                         roundz(as.numeric(UpperCL),digits = 3),")"))  %>%
  mutate(outcome=gsub("BLUE","500mg",outcome),
         outcome=gsub("PINK","1000mg",outcome))


birth_mode <- read.csv(file="../Results/secondary/Mode of delivery 09MAR22.csv",  
                       stringsAsFactors = F) %>%
  rename(outcome = label) %>%
  mutate(Estimate=paste0(roundz(as.numeric(OddsRatioEst),digits = 3), " (", 
                         roundz(as.numeric(LowerCL),digits = 3),", ",
                         roundz(as.numeric(UpperCL),digits = 3),")"))  %>%
  mutate(outcome=gsub("BLUE","500mg",outcome),
         outcome=gsub("PINK","1000mg",outcome))



GA <- obstab(birth_cont, "Delivery_GA Gestational age at delivery  500mg-1000mg", 
             iron, Delivery_GA, cat=0)
GDM <- obstab(birth_or, "GDM gestational diabetes 500mg/1000mg", 
              iron, GDM)
PE <- obstab(birth_or, "PE Pre-eclampsia 500mg/1000mg", 
              iron, PE)
APH <- obstab(birth_or, "APH Antepartum haemorrhage 500mg/1000mg", 
             iron, APH)
PPH <- obstab(birth_or, "PPH Postpartum haemorrhage 500mg/1000mg", 
              iron, PPH)
preterm <- obstab(birth_or, "preterm Preterm Delivery 500mg/1000mg", 
              iron, preterm)

mode <- ordtab(birth_mode, "Mode_Delivery Mode of Delivery 500mg/1000mg", 
                     iron, Delivery) %>%
  mutate(outcome = ifelse(Response == "Instrumental vaginal", outcome, NA_character_),
         Estimate = ifelse(Response == "Normal vaginal", "1.0", Estimate))

birth_tab <- bind_rows(GA, GDM, PE, APH, PPH, preterm, mode)


## infant outcomes

infant_out <- read.csv(file="../Results/secondary/Infant birth outcomes gee 15FEB22.csv",  
                   stringsAsFactors = F) %>%
  mutate(Estimate=paste0(roundz(as.numeric(MeanEstimate),digits = 3), " (", 
                         roundz(as.numeric(MeanLowerCL),digits = 3),", ",
                         roundz(as.numeric(MeanUpperCL),digits = 3),")")) %>%
  rename(p=ProbChiSq) %>%
  mutate(outcome=gsub("BLUE","500mg",outcome),
         outcome=gsub("PINK","1000mg",outcome))

infant_ASQ <- read.csv(file="../Results/secondary/Infant ASQ gee 15FEB22.csv",  
                   stringsAsFactors = F) %>%
  mutate(Estimate=paste0(roundz(as.numeric(MeanEstimate),digits = 3), " (", 
                         roundz(as.numeric(MeanLowerCL),digits = 3),", ",
                         roundz(as.numeric(MeanUpperCL),digits = 3),")")) %>%
  rename(p=ProbChiSq) %>%
  mutate(outcome=gsub("BLUE","500mg",outcome),
         outcome=gsub("PINK","1000mg",outcome))



BWT <- obstab(infant_out, "BWT Birthweight 500mg-1000mg", 
             iron, BWT, cat=0)
BWC <- obstab(infant_out, "BWC Birthweight Centile 500mg-1000mg", 
              iron, BWC, cat=0)
HC <- obstab(infant_out, "HC Head Circumference 500mg-1000mg", 
              iron, HC, cat=0)
Length <- obstab(infant_out, "Length Length 500mg-1000mg", 
             iron, Length, cat=0)

ASQ_COMM_TOTAL <- obstab(infant_ASQ, "ASQ_COMM_TOTAL Communication 500mg-1000mg", 
                 infant, ASQ_COMM_TOTAL, cat=0)
ASQ_GM_TOTAL <- obstab(infant_ASQ, "ASQ_GM_TOTAL Gross Motor 500mg-1000mg", 
                         infant, ASQ_GM_TOTAL, cat=0)
ASQ_FM_TOTAL <- obstab(infant_ASQ, "ASQ_FM_TOTAL Fine motor 500mg-1000mg", 
                       infant, ASQ_FM_TOTAL, cat=0)
ASQ_PS_TOTAL <- obstab(infant_ASQ, "ASQ_PS_TOTAL Problem Solving 500mg-1000mg", 
                       infant, ASQ_PS_TOTAL, cat=0)
ASQ_PERS_TOTAL <- obstab(infant_ASQ, "ASQ_PERS_TOTAL Personal-Social 500mg-1000mg", 
                       infant, ASQ_PERS_TOTAL, cat=0)

infant_tab <- bind_rows(BWT, BWC, HC, Length,
                        ASQ_COMM_TOTAL, ASQ_GM_TOTAL, ASQ_FM_TOTAL,
                        ASQ_PS_TOTAL, ASQ_PERS_TOTAL)
save(iron_sectab, birth_tab, infant_tab, file="../Data/sec_results.Rdata")


#
## MH outcomes
