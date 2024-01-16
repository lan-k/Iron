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
library(car)
library(emmeans)
library(jomo)
library(mitml)
library(mitools)

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


#---- lmer_imp ----

##using jomo package



mh_imp <- function(outvar, basevars, adjvars = basevars, n.imp, n.burn, n.between) {
  
  set.seed=2024
  dat <- mh %>%
    select(all_of(outvar), all_of(basevars),all_of(adjvars),
           rand, stratification, Study_ID, time)
  
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
  
  
  imp.mean3<-list(rep(NA,n.imp))  #6 week estimate
  imp.se3<-list(rep(NA,n.imp))
  imp.mean6<-list(rep(NA,n.imp))  #12 month estimate
  imp.se6<-list(rep(NA,n.imp))
  
  form = formula(paste0(outvar, " ~ rand*time + stratification +", 
                        paste0(adjvars, collapse = "+"), "+ (1|clus)"))
  
  for (k in 1:n.imp) {
    dat.temp<-imp[imp$Imputation==k,]
    
    dat.temp<-data.frame(dat.temp, dat %>%  
                             select(Study_ID, all_of(adjvars)))

    fit <- lmer(formula = form,
                data = dat.temp)
    rt <- emmeans(fit, ~ rand * time )  #estimated means
    (pairs <- summary(pairs(rt, simple="rand")))  #contrasts
    
    
    imp.mean3[[k]]<-as.numeric(pairs[1,"estimate"])  #mean on imputation k 
    imp.se3[[k]]<-as.numeric(pairs[1,"SE"]) #se of imputation k
    imp.mean6[[k]]<-as.numeric(pairs[2,"estimate"])  #mean on imputation k 
    imp.se6[[k]]<-as.numeric(pairs[2,"SE"]) #se of imputation k
  }
  
  #6 weeks
  MIcomb3 <- mitools::MIcombine(imp.mean3,imp.se3)  #combines using Rubin's rules
  estimate3 <- MIcomb3$coefficients
  se3 <- sqrt(MIcomb3$variance)
  p3 <- 2*pnorm(-abs(estimate3/se3))
  
  #12 months
  # MIcomb6 <- miceadds::pool_mi(qhat=list(imp.mean6),se=list(imp.se6))  #combines using Rubin's rules
  
  MIcomb6 <- mitools::MIcombine(imp.mean6,imp.se6)  #combines using Rubin's rules
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
# outvar = "MH_EPDS"
# basevars = c("BL_MH_EPDS", "MH_EPDS_12")
n.burn=as.integer(1000)
n.between=as.integer(1000)
n.imp=as.integer(5)



# HAMD_DEP,HAMA_ANX, SIGHAD, GAF, FAST, EPDS 


(hama_imp <- mh_imp(outvar = "MH_HAMA_ANX", 
                    basevars = c("BL_MH_HAMA_ANX", "Age","Education_years"), 
                    n.imp = n.imp, n.burn = n.burn, n.between = n.between)) 
#singular with income

(hamd_imp <- mh_imp(outvar = "MH_HAMD_DEP", 
                    basevars = c("BL_MH_HAMD_DEP", "Age", "Education_years"), 
                    n.imp = n.imp, n.burn = n.burn, n.between = n.between))
# singular with "income",



(SIGHAD_imp <- mh_imp(outvar = "MH_SIGHAD", 
                    basevars = c("BL_MH_SIGHAD", "Age",  "Education_years"), 
                    n.imp = n.imp, n.burn = n.burn, n.between = n.between))
# singular with , "income"


(GAF_imp <- mh_imp(outvar = "MH_GAF", 
                      basevars = c("BL_MH_GAF", "Age", "income", "Education_years"), 
                   vars = c("BL_MH_GAF", "Age", "income", "Education_years"), 
                      n.imp = n.imp, n.burn = n.burn, n.between = n.between))

(FAST_imp <- mh_imp(outvar = "MH_FAST", 
                   basevars = c("BL_MH_FAST"), 
                   adjvars = c("BL_MH_FAST", "Age", "income", "Education_years"), 
                   n.imp = n.imp, n.burn = n.burn, n.between = n.between))
#singular with "income",, "Age", "Education_years and baseline  -NOT WORKING

(epds_imp <- mh_imp(outvar = "MH_EPDS", 
                    basevars = c("BL_MH_EPDS",  "MH_EPDS_12"), 
                    adjvars = c("BL_MH_EPDS","Age" ,"income" ,  "Education_years"),
                    n.imp = 5, n.burn = n.burn, n.between = n.between))  
#singular ,"income" ,  "Education_years" , "Age"  - NOT WORKING



save(hama_imp, hamd_imp, SIGHAD_imp, GAF_imp, FAST_imp, epds_imp, file="../mh_imp.rds")


#----- specify levels of variables (only relevent for variables
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