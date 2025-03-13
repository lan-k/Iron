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
