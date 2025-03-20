##mediate examples

library(mediation)

# Examples with JOBS II Field Experiment

# **For illustration purposes a small number of simulations are used**

data(jobs)

####################################################
# Example 1: Linear Outcome and Mediator Models
####################################################
b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
c <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs)

# Estimation via quasi-Bayesian approximation
contcont <- mediate(b, c, sims=50, treat="treat", mediator="job_seek")
summary(contcont)
plot(contcont)

# }
# NOT RUN {
# Estimation via nonparametric bootstrap
contcont.boot <- mediate(b, c, boot=TRUE, sims=50, treat="treat", mediator="job_seek")
summary(contcont.boot)
# }
# NOT RUN {
# Allowing treatment-mediator interaction
d <- lm(depress2 ~ treat + job_seek + treat:job_seek + econ_hard + sex + age, data=jobs)
contcont.int <- mediate(b, d, sims=50, treat="treat", mediator="job_seek")
summary(contcont.int)

# Allowing ``moderated mediation'' with respect to age
b.int <- lm(job_seek ~ treat*age + econ_hard + sex, data=jobs)
d.int <- lm(depress2 ~ treat*job_seek*age + econ_hard + sex, data=jobs)
contcont.age20 <- mediate(b.int, d.int, sims=50, treat="treat", mediator="job_seek",
                          covariates = list(age = 20))
contcont.age70 <- mediate(b.int, d.int, sims=50, treat="treat", mediator="job_seek",
                          covariates = list(age = 70))
summary(contcont.age20)
summary(contcont.age70)

# Continuous treatment
jobs$treat_cont <- jobs$treat + rnorm(nrow(jobs))  # (hypothetical) continuous treatment
b.contT <- lm(job_seek ~ treat_cont + econ_hard + sex + age, data=jobs)
c.contT <- lm(depress2 ~ treat_cont + job_seek + econ_hard + sex + age, data=jobs)
contcont.cont <- mediate(b.contT, c.contT, sims=50, 
                         treat="treat_cont", mediator="job_seek",
                         treat.value = 4, control.value = -2)
summary(contcont.cont)

# Categorical treatment 
# }
# NOT RUN {
b <- lm(job_seek ~ educ + sex, data=jobs)
c <- lm(depress2 ~ educ + job_seek + sex, data=jobs)

# compare two categories of educ --- gradwk and somcol
model.cat <- mediate(b, c, treat="educ", mediator="job_seek", sims=50, 
                     control.value = "gradwk", treat.value = "somcol")
summary(model.cat)
# }
# NOT RUN {
######################################################
# Example 2: Binary Outcome and Ordered Mediator
######################################################
# }
# NOT RUN {
jobs$job_disc <- as.factor(jobs$job_disc)
b.ord <- polr(job_disc ~ treat + econ_hard + sex + age, data=jobs,
              method="probit", Hess=TRUE)
d.bin <- glm(work1 ~ treat + job_disc + econ_hard + sex + age, data=jobs,
             family=binomial(link="probit"))
ordbin <- mediate(b.ord, d.bin, sims=50, treat="treat", mediator="job_disc")
summary(ordbin)

# Using heteroskedasticity-consistent standard errors
ordbin.rb <- mediate(b.ord, d.bin, sims=50, treat="treat", mediator="job_disc",
                     robustSE=TRUE)
summary(ordbin.rb)

# Using non-parametric bootstrap
ordbin.boot <- mediate(b.ord, d.bin, sims=50, treat="treat", mediator="job_disc",
                       boot=TRUE)
summary(ordbin.boot)
# }
# NOT RUN {
######################################################
# Example 3: Quantile Causal Mediation Effect
######################################################
require(quantreg)
c.quan <- rq(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs,
             tau = 0.5)  # median
contquan <- mediate(b, c.quan, sims=50, treat="treat", mediator="job_seek")
summary(contquan)

######################################################
# Example 4: GAM Outcome
######################################################
# }
# NOT RUN {
require(mgcv)
c.gam <- gam(depress2 ~ treat + s(job_seek, bs="cr") + 
               econ_hard + sex + age, data=jobs)
contgam <- mediate(b, c.gam, sims=10, treat="treat", 
                   mediator="job_seek", boot=TRUE)
summary(contgam)

# With interaction
d.gam <- gam(depress2 ~ treat + s(job_seek, by = treat) + 
               s(job_seek, by = control) + econ_hard + sex + age, data=jobs)
contgam.int <- mediate(b, d.gam, sims=10, treat="treat", mediator="job_seek",
                       control = "control", boot=TRUE)
summary(contgam.int)
# }
# NOT RUN {
######################################################
# Example 5: Multilevel Outcome and Mediator Models
######################################################
# }
# NOT RUN {
require(lme4)

# educ: mediator group
# occp: outcome group

# Varying intercept for mediator 
model.m <- glmer(job_dich ~ treat + econ_hard + (1 | educ), 
                 family = binomial(link = "probit"), data = jobs)

# Varying intercept and slope for outcome
model.y <- glmer(work1 ~ treat + job_dich + econ_hard + (1 + treat | occp), 
                 family = binomial(link = "probit"), data = jobs)

# Output based on mediator group ("educ")
multilevel <- mediate(model.m, model.y, treat = "treat", 
                      mediator = "job_dich", sims=50, group.out="educ")

# Output based on outcome group ("occp")
multilevel <- mediate(model.m, model.y, treat = "treat", 
                      mediator = "job_dich", sims=50) 

# Group-average effects  
summary(multilevel)

# Group-specific effects organized by effect
summary(multilevel, output="byeffect")
# plot(multilevel, group.plots=TRUE)
# See summary.mediate.mer and plot.mediate.mer for detailed explanations 

# Group-specific effects organized by group
summary(multilevel, output="bygroup")
# See summary.mediate.mer for detailed explanations 
# }



######################################################
# Example 6: Multiple Outcome and Mediator Models
######################################################
# Not run: 
# Hypothetical example

datasets <- list(T1 = T1, T2 = T2)
# List of data frames corresponding to the two different treatment variables 
#"T1vsCont" and "T2vsCont". 
# Each data set has its respective treatment variable.

mediators <- c("M1", "M2") 
# Vector of mediator names, all included in each data frame.

outcome <- c("Ycont1","Ycont2")
# Vector of outcome variable names, again all included in each data frame.

treatment <- c("T1vsCont", "T2vsCont")
# Vector of treatment variables names; must begin with identical strings with dataset 
# names in 'datasets'.

covariates <- c("X1 + X2")
# Set of covariates (in each data set), entered using the standard model formula format.

x <- mediations(datasets, treatment, mediators, outcome, covariates,
                families=c("gaussian","gaussian"), interaction=FALSE, 
                conf.level=.90, sims=50) 
# Runs 'mediate' iteratively for each variable combinations, with 'lm' on both mediator 
# and outcome model.

summary(x)  # tabular summary of results for all model combinations
plot(x)  # graphical summary of results for all model combinations at once

plot(x$Ycont1.T1vsCont.M1) 
# Individual 'mediate' outputs are stored as list elements and can be 
# accessed using the usual "$" operator.

## End(Not run)


######################################################
# https://rpubs.com/mbounthavong/mediation_analysis_r
######################################################


### step 3: Load the MEPS package
library("MEPS") ## You need to load the library every time you restart R

### Step 4: Load the other libraries
library("dplyr")          # Data wrangling
library("gtsummary")      # Create tables
library("expss")          # Relabel variables
library("bda")            # Perform Sobel text
library("mediation")        # Perform the bootstrap approach

# There are two ways to load data from AHRQ MEPS website:
#### Method 1: Load data from AHRQ MEPS website
hc2021 = read_MEPS(file = "h233")

#### Method 2: Load data from AHRQ MEPS website
hc2021 = read_MEPS(year = 2021, type = "FYC")

### Step 5: Change column names to lowercase
names(hc2021) <- tolower(names(hc2021))

### Step 6: Select specific variables
### 2021
hc2021p = hc2021 %>%
  rename(
    workdays = ddnwrk21,
    diabetes = diabdx_m18,
    health_status = rthlth31) %>%
  dplyr::select(                        ## NOTE: there is a weird issue when you don't involve dplyr::select bc MASS is preferred (URL: https://stackoverflow.com/questions/48161431/select-statement-error-unused-argument)
    dupersid, 
    workdays, 
    diabetes, 
    health_status,
    sex)
hc2021p$year <- 2021

### Step 7: Clean data (We don't want to include any missing or NA responses)
hc2021p = hc2021p %>%
  filter(workdays >= 0,
         diabetes >= 1,
         health_status >= 1)

# We want "No diabetes" to have a value of 0 because it will make interpreting the model easier
hc2021p$diabetes[hc2021p$diabetes == 2] = 0

### Step 8: Convert to factor and add labels
hc2021p$sex <- factor(hc2021p$sex,
                      levels = c(1, 2),
                      labels = c("Male", "Female"))

hc2021p$health_status <- factor(hc2021p$health_status, 
                                levels = c(1, 2, 3, 4, 5),
                                labels = c("Excellent", "Very good", "Good", "Fair", "Poor"))

hc2021p$diabetes <- factor(hc2021p$diabetes,
                           levels = c(0, 1),
                           labels = c("No diabetes", "Diabetes"))

### Step 9: Relabel the variables
hc2021p <- apply_labels(hc2021p,
                        diabetes = "Diabetes status", 
                        workdays = "Days missed from work", 
                        health_status = "Perceived health status", 
                        sex = "Sex")



hc2021p %>%
  tbl_summary(include = c(health_status, sex, workdays), 
              by = diabetes, 
              statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  modify_header(label = "**Variable***") %>%
  bold_labels()


direct.model <- glm(workdays ~ diabetes, data = hc2021p, family = gaussian(link = "identity"))
round(cbind(coef(direct.model), confint(direct.model)), 3)

tbl_regression(direct.model, estimate_fun = ~ style_number(.x, digits = 3))


indirect.model1 <- glm(as.numeric(health_status) ~ diabetes, data = hc2021p, family = gaussian(link = "identity")) ## We need to convert the `health_status` variable to numeric to make sense of the linear regression model.
round(cbind(coef(indirect.model1), confint(indirect.model1)), 3)

tbl_regression(indirect.model1, estimate_fun = ~ style_number(.x, digits = 3))

indirect.model2 <- glm(workdays ~ diabetes + as.numeric(health_status), data = hc2021p, family = gaussian(link = "identity"))
round(cbind(coef(indirect.model2), confint(indirect.model2)), 3)

tbl_regression(indirect.model2, estimate_fun = ~ style_number(.x, digits = 3))


t1 <- tbl_regression(direct.model, estimate_fun = ~ style_number(.x, digits = 3)) %>% modify_column_hide(columns = p.value) %>% modify_column_hide(ci)
t2 <- tbl_regression(indirect.model1, estimate_fun = ~ style_number(.x, digits = 3)) %>% modify_column_hide(columns = p.value) %>% modify_column_hide(ci)
t3 <- tbl_regression(indirect.model2, estimate_fun = ~ style_number(.x, digits = 3)) %>% modify_column_hide(columns = p.value) %>% modify_column_hide(ci)

tbl_merge(
  tbls = list(t1, t2, t3),
  tab_spanner = c("**Total effect**", "**First-part indirect effect**", "**Second-part direct + indirect effects**")
)


## mediation.test(mv,iv,dv), where mv is the mediator variable, 
# iv is the independent varible, and dv is the dependent variable)
mediation <- mediation.test(as.numeric(hc2021p$health_status), hc2021p$diabetes, hc2021p$workdays)


## Note: We have to convert health_status to a numeric since errors will occur when using it as a factor.
hc2021p$health_status <- as.numeric(hc2021p$health_status)

## First-part indirect effect model:
indirect.model1 <- glm(health_status ~ diabetes, data = hc2021p, family = gaussian(link = "identity"))

## Second-part direct + indirect effect model:
indirect.model2 <- glm(workdays ~ diabetes + health_status, data = hc2021p, family = gaussian(link = "identity"))

## Mediation analysis with 1000 simulations
mediation.results <- mediate(indirect.model1, indirect.model2, treat = 'diabetes', mediator = 'health_status', boot = TRUE, sims = 1000)
summary(mediation.results)
round(mediation, 3)