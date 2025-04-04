##functions
roundz = function(x, digits){
  dformat = paste('%.', digits, 'f', sep='')
  x = sprintf(dformat, round(x, digits))
  return(x)
}


quantile_str <- function(x, probs = c(0.25, 0.5, 0.75), digits=0) {
  value = quantile(x, probs, na.rm=T, digits=digits)
  n=sum(!is.na(x))
  
  med_iqr = paste0(as.character(roundz(value[2],digits=digits)), " (", 
                   as.character(roundz(value[1],digits=digits)), "-",
                   as.character(roundz(value[3],digits=digits)), ")",
                   " (N=",n,")")
  
  return(med_iqr)
}


mean_sd <- function(x,  digits=2, exp = F) {
  
  if (exp) {
    mean = exp(mean(x, na.rm=T))
    sd = exp(sd(x, na.rm=T))
    n=sum(!is.na(x))
    
  } else {
    mean = mean(x, na.rm=T)
    sd = sd(x, na.rm=T)
    n=sum(!is.na(x))
    
  }
  
  
  mean_sd = paste0(as.character(roundz(mean,digits=digits)), " (", 
                   as.character(roundz(sd,digits=digits)), ")",
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

est_ci = function(est, ci, digits=2) {
  return(paste0(roundz(est, digits), " (",
                roundz(ci[1], digits), ", ",
                roundz(ci[2], digits), ")"))
}

fit_censreg = function(form, df, leftlim = -Inf, rightlim = Inf,
                       label = "") {
  
  df = df %>% mutate(time = factor(time))
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
  
  df = df %>% mutate(time = factor(time))
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
                       df = mh_il6_vitd, timept = NULL,exp=F,
                        nsim = 500) {
  dat <- df %>%
    dplyr::select(study_id, treat, time, stratification, age, income,
                  education_years,
                   all_of(c(outcome, outcome_base, 
                  mediator, mediator_base ))) %>%
    na.omit() %>%
    arrange(study_id, time) %>%
    mutate(timept=ifelse(as.character(time == "1"), "6 weeks", "12 months"),
           time=factor(time))
  
  ##descriptives
  if (exp) {
    
    desc = dat %>%
      group_by(treat, timept) %>%
      reframe(across(all_of(c(mediator, outcome)), \(x) mean_sd(x, digits=2, exp=T))) %>% 
      ungroup()
    
    
    base = dat %>%
      dplyr::select(treat, all_of(c(mediator_base, outcome_base))) %>%
      mutate(timept = "Screening") %>%
      group_by(treat, timept) %>%
      reframe(across(all_of(c( mediator_base)), \(x) mean_sd(x, digits=2, exp=T)),
              across(all_of(c( outcome_base)), \(x) mean_sd(x, digits=2, exp=F))) %>% 
      ungroup() 
    
    
  } else {
    
    desc = dat %>%
      group_by(treat, timept) %>%
      reframe(across(all_of(c(mediator, outcome)), \(x) mean_sd(x, digits=2))) %>% 
      ungroup()
    
    
    base = dat %>%
      dplyr::select(treat, all_of(c(mediator_base, outcome_base))) %>%
      mutate(timept = "Screening") %>%
      group_by(treat, timept) %>%
      reframe(across(all_of(c( mediator_base, outcome_base)), \(x) mean_sd(x, digits=2))) %>% 
      ungroup() 
    
    
  }

  
  
  colnames(base) = c("treat","timept", mediator, outcome)
  
  
  base_wide = base %>%
    pivot_wider(names_from = "treat", 
                values_from = c(mediator, outcome) )

  
  desc_wide <- desc %>%
    dplyr::select(treat,timept, all_of(c(mediator, outcome))) %>%
    pivot_wider(names_from = "treat", 
                values_from = c(mediator, outcome) ) %>% 
    arrange(desc(timept))
  
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
  if (!is.null(timept)) {
    contcont <- mediate(fit_med, outfit, sims=nsim, treat="treat", 
                        mediator=mediator, group.out = "study_id",  
                        covariates = list(time = timept))
  } else {
    contcont <- mediate(fit_med, outfit, sims=nsim, treat="treat", 
                        mediator=mediator, group.out = "study_id")
  }
  
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
                  acme0=est_ci(s$d0, s$d0.ci),
                  acme0_p = format.pval(s$d0.p, eps=0.001, digits=2),
                  acme1=est_ci(s$d1,s$d1.ci), 
                  acme1_p = format.pval(s$d1.p, eps=0.001, digits=2),
                  acme = est_ci(s$d.avg,s$d.avg.ci),
                  acme_p = format.pval(s$d.avg.p, eps=0.001, digits=2),
                  ade0=est_ci(s$z0,s$z0.ci),
                  ade0_p = format.pval(s$z0.p, eps=0.001, digits=2),
                  ade1=est_ci(s$z1,s$z1.ci),
                  ade1_p = format.pval(s$z1.p, eps=0.001, digits=2),
                  ade=est_ci(s$z.avg,s$z.avg.ci),
                  ade_p = format.pval(s$z.avg.p, eps=0.001, digits=2),
                  ate=est_ci(s$tau.coef,s$tau.ci),
                  ate_p = format.pval(s$tau.p, eps=0.001, digits=2),
                  prop_med0 = est_ci(s$n0,s$n0.ci,digits=3),
                  prop_med1 = est_ci(s$n1,s$n1.ci,digits=3),
                  prop_med = est_ci(s$n.avg,s$n.avg.ci,digits=3)
                  )
  
  rownames(df) = c()
  
  return(df)
  
}

