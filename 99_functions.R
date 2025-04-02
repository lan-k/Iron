##functions

quantile_str <- function(x, probs = c(0.25, 0.5, 0.75), digits=0) {
  value = quantile(x, probs, na.rm=T, digits=digits)
  n=sum(!is.na(x))
  
  med_iqr = paste0(as.character(round(value[2],digits=digits)), " (", 
                   as.character(round(value[1],digits=digits)), "-",
                   as.character(round(value[3],digits=digits)), ")",
                   " (N=",n,")")
  
  return(med_iqr)
}


mean_sd <- function(x,  digits=2) {
  
  mean = mean(x, na.rm=T)
  sd = sd(x, na.rm=T)
  n=sum(!is.na(x))
  
  mean_sd = paste0(as.character(round(mean,digits=digits)), " (", 
                   as.character(round(sd,digits=digits)), ")",
                   " (N=",n,")")
  
  return(mean_sd)
}

mean_log <- function(x,  digits=2) {
  
  mean = exp(mean(log(x), na.rm=T))
  sd = exp(sd(log(x), na.rm=T))
  n=sum(!is.na(x))
  
  mean_sd = paste0(as.character(round(mean,digits=digits)), " (", 
                   as.character(round(sd,digits=digits)), ")",
                   " (N=",n,")")
  
  return(mean_sd)
}



fit_censreg = function(form, df, leftlim = -Inf, rightlim = Inf,
                       label = "") {
  
  
  fit <- censReg( form,  data = df, method = "BHHH",
                  left = leftlim, right=rightlim)
                             
  summary( fit )
  
  coef <- coef(fit)

  #drop last 2 columns
  ncol=length(coef)-2
  coef <- coef(fit)[1:ncol]
  vcov <- vcov(fit)[1:ncol, 1:ncol]
  
  #create a grid because emmeans doesn't work with censReg
  rg <- emmeans::qdrg(form, 
                      data=df, coef=coef, vcov=vcov, 
                      df=fit$df.residual)
  
  # emmeans(rg, pairwise ~ rand|time, adjust="none")
  p <- summary(pairs(emmeans(rg, ~ rand * time ), 
                      simple="rand"))  #contrasts PINK - BLUE
  
  
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
  
  if (grepl("log", as.character(form)[2])) {
    outdf = data.frame(Outcome = c(label,""),
                       Contrast = c("500mg vs 1000mg",""), 
                       Level=time, Estimate = exp(-estimate), 
                       lower = exp(-estimate - 1.96*se),
                       upper = exp(-estimate + 1.96*se),
                       se = se, p=p)
    
  } else {
    outdf = data.frame(Outcome = c(label,""),
                       Contrast = c("500mg - 1000mg",""), 
                       Level=time, Estimate = -estimate, 
                       lower = -estimate - 1.96*se,
                       upper = -estimate + 1.96*se,
                       se = se, p=p)
  }
  
  
  
  rownames(outdf) <- c()
  
  return(outdf)
  
}

fit_lmer = function(form, df, label = "") {
  
  fit <- lmer(form , data = df) 
  

  # anova(fit)
  p <- summary(pairs(emmeans(fit, ~ rand * time ), 
                      simple="rand"))  #contrasts PINK - BLUE
  
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
  
  if (grepl("log", as.character(form)[2])) {
    outdf = data.frame(Outcome = c(label,""),
                        Contrast = c("500mg vs 1000mg",""), 
                        Level=time, Estimate = exp(-1*estimate), 
                        lower = exp(-estimate - 1.96*se),
                        upper = exp(-estimate + 1.96*se),
                        se = se, p=p)
    
  } else {
    outdf = data.frame(Outcome = c(label,""),
                        Contrast = c("500mg-1000mg",""), 
                        Level=time, Estimate = -1*estimate, 
                        lower = -estimate - 1.96*se,
                        upper = -estimate + 1.96*se,
                        se = se, p=p)
    
  }
  
  
  
  rownames(outdf) <- c()
  
  return(outdf)
  
}

lmer_nolim = function(outcome, outlab, screen, df, limit) {
  
  df = df %>%
    mutate(out={{outcome}},
           out =ifelse(out == limit & limit > 0, limit*runif(1), out),
           screen={{screen}},
           screen =ifelse(screen == limit & limit > 0, limit*runif(1), screen))
  
  form = formula(log(out) ~ rand*time + stratification + log(screen) +
                   (1|study_id))

  fit <- lmer(formula = form, data = df) 
  summary(fit) 
  
  (p <- summary(pairs(emmeans(fit, ~ rand * time ), 
                      simple="rand")))  #contrasts PINK - BLUE i.e. 1000 - 500
  
  
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
  
  out_df = data.frame(Outcome = c(outlab,""),
                      Contrast = c("500mg vs 1000mg",""), 
                      Level=time, Estimate = exp(-estimate), 
                      lower = exp(-estimate - 1.96*se),
                      upper = exp(-estimate + 1.96*se),
                      se = se, p=p)
  
  rownames(out_df) <- c()
  
  return(out_df)
  
}

mh_mediate = function( outcome, outcome_base, mediator, mediator_base,
                       df = mh_il6_vitd,
                        nsim = 500) {
  dat <- df %>%
    dplyr::select(study_id, treat, time, stratification, age, income,
                  education_years,
                   all_of(c(outcome, outcome_base, 
                  mediator, mediator_base ))) %>%
    na.omit() %>%
    mutate(timept=ifelse(as.character(time == "1"), "FU1", "FU2"))
  
  ##descriptives
  

  
  desc = dat %>%
    group_by(treat, timept) %>%
    reframe(across(all_of(c(mediator, outcome)), \(x) mean_sd(x, digits=2))) %>% 
    ungroup()
  

  base = dat %>%
    dplyr::select(treat, all_of(c(mediator_base, outcome_base))) %>%
    mutate(timept = "SCR") %>%
    group_by(treat, timept) %>%
    reframe(across(all_of(c( mediator_base, outcome_base)), \(x) mean_sd(x, digits=2))) %>% 
    ungroup() 
  
  
  colnames(base) = c("treat","timept", mediator, outcome)
  
  
  base_wide = base %>%
    pivot_wider(names_from = "treat", 
                values_from = c(mediator, outcome) )

  
  desc_wide <- desc %>%
    dplyr::select(treat,timept, all_of(c(mediator, outcome))) %>%
    pivot_wider(names_from = "treat", 
                values_from = c(mediator, outcome) )
  
  desc_wide = bind_rows(base_wide, desc_wide)
  
  
  #mediation
  
  base_form = "treat* time + stratification + age + income + education_years + (1|study_id)"
  mform = formula(paste0(mediator,"~", mediator_base, "+", base_form))
  
  
  fit_med <- lme4::lmer(formula = mform, data = dat)
  
  # outcome
  outform = formula(paste0(outcome, "~ ", base_form, "+",
                    paste(c(outcome_base, mediator, mediator_base), collapse="+") ) )
  
  
  outfit <- lme4::lmer(formula = outform, data = dat)
  
  
  ##mediation
  contcont <- mediate(fit_med, outfit, sims=nsim, treat="treat", 
                      mediator=mediator, group.out = "study_id")
  return(list(desc = desc_wide, medfit = contcont))
}


mh_mediate_lm = function(outcome, outcome_base, mediator, mediator_base,
                         df= mh_il6_vitd , timept = 1,
                      nsim = 500) {
  dat <- df %>%
    filter(time==timept) %>%
    dplyr::select(study_id, time, treat, stratification, age, income,
                  education_years,
                  all_of(c(outcome, outcome_base, 
                           mediator, mediator_base ))) %>%
    na.omit() 
  
  
  base_form = "treat + stratification + age + income + education_years "
  mform = formula(paste0(mediator,"~", mediator_base, "+", base_form))
  
  
  fit_med <- lm(formula = mform, data = dat)
  
  # outcome
  outform = formula(paste0(outcome, "~ ", base_form, "+",
                           paste(c(outcome_base, mediator, mediator_base), collapse="+") ) )
  
  
  outfit <- lm(formula = outform, data = dat)
  
  
  ##mediation
  contcont <- mediate(fit_med, outfit, sims=nsim, treat="treat", 
                      mediator=mediator)
  return(contcont)
}


tab_mh_mediate = function(medfit, Outcome = NULL, mediator = NULL) {
  
 
  s = summary(medfit)
  
  df = data.frame(Outcome = Outcome,
                  Mediator = mediator,
                  acme0=s$d0,
                  acme0_lo = s$d0.ci[1],
                  acme0_hi = s$d0.ci[2],
                  acme0_p = s$d0.p,
                  acme1=s$d1, 
                  acme1_lo = s$d1.ci[1],
                  acme1_hi = s$d1.ci[2],
                  acme1_p = s$d1.p,
                  acme = s$d.avg,
                  acme_lo = s$d.avg.ci[1],
                  acme_hi = s$d.avg.ci[2],
                  acme_p = s$d.avg.p,
                  ade0=s$z0,
                  ade0_lo = s$z0.ci[1],
                  ade0_hi = s$z0.ci[2],
                  ade0_p = s$z0.p,
                  ade1=s$z1,
                  ade1_lo = s$z1.ci[1],
                  ade1_hi = s$z1.ci[2],
                  ade1_p = s$z1.p,
                  ade=s$z.avg,
                  ade_lo = s$z.avg.ci[1],
                  ade_hi = s$z.avg.ci[2],
                  ade_p = s$z.avg.p,
                  ate=s$tau.coef,
                  ate_lo = s$tau.ci[1],
                  ate_hi = s$tau.ci[2],
                  ate_p = s$tau.p,
                  prop_med0 = s$n0,
                  prop_med0_lo = s$n0.ci[1],
                  prop_med0_hi = s$n0.ci[2],
                  prop_med1 = s$n1,
                  prop_med1_lo = s$n1.ci[1],
                  prop_med1_hi = s$n1.ci[2],
                  prop_med = s$n.avg,
                  prop_med_lo = s$n.avg[1],
                  prop_med_hi = s$n.avg[2]
                  )
  
  rownames(df) = c()
  
  return(df)
  
}
