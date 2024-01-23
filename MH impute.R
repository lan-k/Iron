#impute MH data using https://www.tandfonline.com/doi/full/10.1080/10543406.2013.834911
# https://search.r-project.org/CRAN/refmans/mlmi/html/refBasedCts.html
# use ml.lmer function in the R package ‘miceadds’
# HAMD_DEP,HAMA_ANX, SIGHAD, GAF, FAST, EPDS 
# change_MINI_coll n_change_MINI

#---- mh_data ----
rm(list=ls())
library(dplyr)
library(stringr)
library(haven)
library(mice)
library(miceadds)
library(lme4)
# library(car)
library(emmeans)
library(jomo)
# library(mitml)
library(mitools)



n.burn=as.integer(1000)
n.between=as.integer(1000)
n.imp=as.integer(20)


mh <- read_sas(data_file = "../iron_mh.sas7bdat")

mh <- mh %>% 
  mutate(time = factor(time),
         rand=relevel(factor(rand), ref="PINK"),
         stratification=factor(stratification),
         income=factor(ifelse(income==6,5,income)))

levels(mh$time)
table(mh$income,mh$rand)

#---- lmer_epds ----

epds <- lmer(formula = MH_EPDS ~ rand*time + stratification + BL_MH_EPDS + (1|Study_ID),
             data = mh)

summary(epds)
# Anova(epds, test="F" , method = "REML")

rt <- emmeans(epds, ~ rand * time )  #estimated means
(pair_epds <- summary(pairs(rt, simple="rand")))  #contrasts
pair_epds[1,"estimate"]; pair_epds[1,"SE"];pair_epds[1,"p.value"]  #6 weeks
pair_epds[2,"estimate"]; pair_epds[2,"SE"];pair_epds[2,"p.value"]  #12 months



#---- imp_all ----

#impute all the variables at once using jomo
mh_jomo <- function(dat, outvar, basevars,  n.imp, n.burn, n.between) {

  set.seed=2024
  
  
  x <- dat %>%   #covariates at  level 1 with no missing data allowed
    select(time)   %>%     
    mutate(cons=1)  %>%  #for intercept
    data.frame()
  
  x2 <- dat %>% #cluster level covariates in model  
    select(rand, stratification) %>%
    mutate(cons=1)  %>%  #for intercept
    data.frame()
  
  
  y <- dat %>% select(all_of(outvar))  %>% data.frame() #time level variables to be imputed
  y2 <- dat %>% select(all_of(basevars)) %>% data.frame()  #cluster level variables to be imputed
  
  # z = dat %>% #cluster level covariates in model  
  #   select(rand) %>% data.frame()
  
  
  clust = dat$Study_ID
  
  set.seed(1569)
  
  imp <- jomo(Y=y,Y2=y2, X=x,X2=x2, #Z=z, 
              clus=clust, nimp=n.imp, nburn=n.burn, 
              nbetween=n.between, meth="random") #output=0,
  


  return(imp)

}





dat <- mh %>% select(!contains("MINI"))
outvars <- colnames(dat %>% select(contains("MH_")))
basevars <- colnames(dat %>% select(contains("BL_"), 
                                    Age, income, Education_years))


mh_imp_all <- mh_jomo(dat = dat, outvar = outvars, basevars = basevars,  
                       n.imp = n.imp, n.burn=n.burn, n.between=n.between)

save(mh_imp_all, file="../mh_impall.rds")

#---- lmer_pool ----

load(file="../mh_impall.rds")

lmer_pool <- function(imp, outvar, adjvars, n.imp) {
  
  imp.mean3<-list(rep(NA,n.imp))  #6 week estimate
  imp.var3<-list(rep(NA,n.imp))
  imp.mean6<-list(rep(NA,n.imp))  #12 month estimate
  imp.var6<-list(rep(NA,n.imp))
  
  form = formula(paste0(outvar, " ~ rand*time + stratification +", 
                        paste0(adjvars, collapse = "+"), "+ (1|clus)")) #
  
  for (k in 1:n.imp) {
    dat.temp<-imp[imp$Imputation==k,]
    
    dat.temp<-data.frame(dat.temp, dat %>%  
                           select(Study_ID, all_of(adjvars)))
    
    dat.temp$Study_ID = factor(dat.temp$Study_ID)
    
    fit <- lmer(formula = form,data = dat.temp)
    
    rt <- emmeans(fit, ~ rand * time )  #estimated means
    (pairs <- summary(pairs(rt, simple="rand")))  #contrasts
    
    
    imp.mean3[[k]]<-as.numeric(pairs[1,"estimate"])  #mean on imputation k at 6 weeks
    imp.var3[[k]]<-as.numeric(pairs[1,"SE"])^2 #var of imputation k at 6 weeks
    imp.mean6[[k]]<-as.numeric(pairs[2,"estimate"])  #mean on imputation k at 12m
    imp.var6[[k]]<-as.numeric(pairs[2,"SE"])^2 #var of imputation k at 12m
  }
  
  #6 weeks
  MIcomb3 <- mitools::MIcombine(imp.mean3,imp.var3)  #combines using Rubin's rules
  estimate3 <- MIcomb3$coefficients
  se3 <- sqrt(MIcomb3$variance)
  p3 <- 2*pnorm(-abs(estimate3/se3))
  
  #12 months
  # MIcomb6 <- miceadds::pool_mi(qhat=list(imp.mean6),se=list(imp.se6))  #combines using Rubin's rules
  
  MIcomb6 <- mitools::MIcombine(imp.mean6,imp.var6)  #combines using Rubin's rules
  estimate6 <- MIcomb6$coefficients
  se6 <- sqrt(MIcomb6$variance)
  p6 <- 2*pnorm(-abs(estimate6/se6))
  
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


(hama_imp <- lmer_pool(imp= mh_imp_all,outvar = "MH_HAMA_ANX", 
                    adjvars = c("BL_MH_HAMA_ANX", "Age", "income", "Education_years"), 
                    n.imp = n.imp)) 
#singular with income

(hamd_imp <- lmer_pool(imp= mh_imp_all,outvar = "MH_HAMD_DEP", 
                    adjvars = c("BL_MH_HAMD_DEP", "Age", "income", "Education_years"), 
                    n.imp = n.imp))
# singular with "income",



(SIGHAD_imp <- lmer_pool(imp= mh_imp_all,outvar = "MH_SIGHAD", 
                      adjvars = c("BL_MH_SIGHAD", "Age", "income", "Education_years"), 
                      n.imp = n.imp))
# singular with , "income"


(GAF_imp <- lmer_pool(imp= mh_imp_all,outvar = "MH_GAF", 
                   adjvars = c("BL_MH_GAF", "Age", "income", "Education_years"), 
                   n.imp = n.imp))

(FAST_imp <- lmer_pool(imp= mh_imp_all,outvar = "MH_FAST", 
                    adjvars = c("BL_MH_FAST", "Age", "income", "Education_years"), 
                    n.imp = n.imp))

(epds_imp <- lmer_pool(imp= mh_imp_all,outvar = "MH_EPDS", 
                    adjvars = c("BL_MH_EPDS","Age" ,"income" ,  "Education_years"),
                    n.imp = n.imp))  


#BL_total_mini is baseline variable for MINI

imp_pool_results <- bind_rows(hama_imp ,
                              hamd_imp ,
                              SIGHAD_imp ,
                              GAF_imp ,
                              FAST_imp ,
                              epds_imp ) %>%
  mutate(lower=Estimate - 1.96*se,
         upper=Estimate + 1.96*se,
         est_ci = paste0(round(-1*Estimate, 2), " (",
                         round(-1*upper, 2), ", ",
                         round(-1*lower, 2), ")"),
         p = format.pval(p, eps=0.001, digits=2))

write.csv(imp_pool_results, file="../MH Results/MH imputed.csv", row.names=F)

save(imp_pool_results, file="../mh_pool.rds")


#---- lmer_pool_tab ----
library(flextable)

load(file="../mh_pool.rds")


imp_pool_results %>%
  select(Outcome, time, est_ci, p) %>%
  flextable %>%
  autofit() %>%
  theme_box() %>%
  set_header_labels(est_ci = "500mg - 1000mg (95% CI)",
                    time= "Time Point")



#----- specify levels of variables (only relevant for variables
#      with missing values)

## time is id1
## Study_ID is id2
## rand is id3 or ignore


dat <- mh %>%
  select(contains("EPDS"), rand, stratification, Study_ID, time) %>%
  # select(!MH_EPDS_12) %>%
  mutate(time=paste0(Study_ID, time))


variables_levels <- miceadds:::mice_imputation_create_type_vector( colnames(dat), value="")
# leave variables at lowest level blank (i.e., "")
variables_levels[ c("BL_MH_EPDS","MH_EPDS_12","stratification") ] <- "Study_ID"


#----- specify predictor matrix
predmat <- mice::make.predictorMatrix(data=dat)
predmat[, c("Study_ID","rand") ] <- 0
# set -2 for cluster identifier for Study ID variables
# because "2lonly" function is used
predmat[ c("BL_MH_EPDS","MH_EPDS_12"), "Study_ID" ] <- -2


#----- specify imputation methods
method <- mice::make.method(data=dat)
method[c("MH_EPDS")] <- "ml.lmer"
method[c("BL_MH_EPDS","MH_EPDS_12")] <- "2lonly.norm"

#----- specify hierarchical structure of imputation models
levels_id <- list()

#** hierarchical structure for variables
levels_id[[("MH_EPDS")]] <- c("Study_ID", "rand")



#----- specify random slopes
random_slopes <- list()
#** random slopes for variable x1
random_slopes[["MH_EPDS"]] <- list( "Study_ID" =c("BL_MH_EPDS","MH_EPDS_12") )
# if no random slopes should be specified, the corresponding entry can be left empty
# and only a random intercept is used in the imputation model

#----- imputation in mice
imp1 <- mice::mice( dat, maxit=10, m=5, method=method,
                    predictorMatrix=predmat, levels_id=levels_id,
                    # random_slopes=random_slopes,
                    variables_levels=variables_levels )
summary(imp1)
a <- complete(imp1)


epds <- lmer(formula = MH_EPDS ~ rand*time + stratification + BL_MH_EPDS + (1|Study_ID),
             data = a)

summary(epds)
Anova(epds, test="F" , method = "REML")

rt <- emmeans(epds, ~ rand * time )  #estimated means
pair_epds <- summary(pairs(rt, simple="rand"))  #contrasts
pair_epds[1,"estimate"]; pair_epds[1,"SE"];pair_epds[1,"p.value"]  #6 weeks
pair_epds[2,"estimate"]; pair_epds[2,"SE"];pair_epds[2,"p.value"]  #12 months