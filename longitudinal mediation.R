## mediation using method in https://pubmed.ncbi.nlm.nih.gov/32458988/
## web material https://academic.oup.com/aje/article-pdf/189/11/1427/34045752/kwaa092.pdf

First read the data into a matrix/dataframe/vector named mydata
mydata<-read.csv(file="mentalhealth.csv",header=TRUE,sep=",")
N<-nrow(mydata)
#creating an Id variable
mydata$id<-1:N
#L0 has five variables hence L01 L02 L03 L04 L05.
# M1 M2 M3 M4 mediator observed at times 0 1 2 and 3.
head(mydata, n=10L)
Data could not be provided due to confidentiality issues.
#We have mediator measured at four time points.
#To implement the method described in the paper
#we duplicate the data set twice for each mediator,
mydata1<-mydata
mydata2<-mydata
mydata1$Astar<-mydata$A
mydata2$Astar<-1-mydata$A

#thus creating 8 (4x2) replicates of original data,
#we name this duplicated data as newmydata;
newmydata<-rbind(mydata1,mydata2,mydata3,mydata4,mydata5,mydata6,mydata7,mydata8)


# Step-1: We first load the data set (e.g. mental health data) into a vector named “mydata” in R. Then fit a
# generalised linear regression to the mediator (mental psychological distress) at each of the four
# time points separately, conditioning on the exposure A and the confounders listed above (which include
#                                                                                          the history of outcomes and mediators at each time). Since the mediator is a categorical variable we used
# the vector generalised linear model (vglm), with family being multinomial, from the VGAM library in
# R.
m1creg <- vglm(M1~A+factor(L01)+factor(L02)+factor(L03)+factor(L04)+factor(L05)
               ,data=mydata,family=multinomial)
m2creg <- vglm(M2~factor(M1)+A+Y1+factor(L01)+factor(L02)+factor(L03)
               +factor(L04)+factor(L05),data=mydata,family=multinomial)
m3creg <- vglm(M3~factor(M2)+factor(M1)+Y2+Y1+A+factor(L01)+factor(L02)
               +factor(L03)+factor(L04)+factor(L05),data=mydata,family=multinomial)
m4creg <- vglm(M4~factor(M3)+factor(M2)+factor(M1)+Y2+Y1+Y3+A+factor(L01)
               +factor(L02)+factor(L03)+factor(L04)+factor(L05),
               data=mydata,family=multinomial)
# Step-2: In this step we create a new data set by replicating the original two times. Next we create a new
# variable Id, which is a subject identifier, and Astar which takes once takes the value contained in A for
# each subject, and once the opposite value 1-A.
N <- nrow(mydata)
mydata$id <- 1:N
mydata1 <- mydata
mydata2 <- mydata

mydata1$Astar <- mydata$A
mydata2$Astar <- 1-mydata$A
newmydata <- rbind(mydata1,mydata2)
# Step-3: The weights are computed from the predicted probabilities of the above multinomial regressions.
# For the numerator computation in weight equation 11 we used 1􀀀A and for the denominator we use used
# A as explanatory variables:
  #denominator weights
m1pdr <- as.matrix(predict(m1reg,type="response",newdata=mydata))[cbind(1:nrow(mydata),mydata$M1)]
m2pdr <- as.matrix(predict(m2reg,type="response",newdata=mydata))[cbind(1:nrow(mydata),mydata$M2)]
m3pdr <- as.matrix(predict(m3reg,type="response",newdata=mydata))[cbind(1:nrow(mydata),mydata$M3)]
m4pdr <- as.matrix(predict(m4reg,type="response",newdata=mydata))[cbind(1:nrow(mydata),mydata$M4)]


#numerator weights
Dattemp <- mydata
Dattemp$A <- 1 - Dattemp$A
m1pnr<-as.matrix(predict(m1reg,type="response",newdata=Dattemp))[cbind(1:nrow(mydata),mydata$M1)]
m2pnr<-as.matrix(predict(m2reg,type="response",newdata=Dattemp))[cbind(1:nrow(mydata),mydata$M2)]
m3pnr<-as.matrix(predict(m3reg,type="response",newdata=Dattemp))[cbind(1:nrow(mydata),mydata$M3)]
m4pnr<-as.matrix(predict(m4reg,type="response",newdata=Dattemp))[cbind(1:nrow(mydata),mydata$M4)]



# Step-4: The weights for the exposure were created using the logistic regression conditioned on the
# baseline confounders.
#Unstabilized weights
Aregdr <- glm(A~factor(L01)+factor(L02)+factor(L03)+factor(L04)+factor(L05)
              ,data=mydata,family=binomial)
ap <- predict(Areg,type="response")
wa <- ifelse(mydata$A==1,1/ap,1/(1-ap)) #weight of A.

# Step-5: Finally the logistic natural effects model for the dichotomous outcome can be fitted. To obtain
# sandwich standard errors that account for the repeated measures nature of the data, we use the geeglm
# function from the package geepack using the Id variable to indicate dependence. Note that the geeglm
# function requires that the observations be sorted by the Id variable.

library(geepack)
newmydata <- newmydata[order(newmydata$id),]
fit_out <- geeglm(Y~A+Astar+t+I(t*A)+I(T*Astar),data=newmydata,corstr="independence",
                  family="binomial", weights=W,id=newmydata$id,scale.fix=T)

summary(fit_out)


# Code for computing the variance
# The below code is for computing the variance using the method described in section-4 of the manuscript.
# Code for drawing samples
# This section has three components; 1) adding error to coecients of exposure regression, 2) adding error
# to the coecients of mediator regression and 3) computing probabilities of the multinomial regression
#Part-1: generating new set of weights for the exposure regression
library(mvtnorm)
awts<-function(regfit){
  sreg<-summary(regfit)
  cv<-sreg$cov.unscaled
  mdl<-model.matrix(regfit)
  cf<-sreg$coefficients[,1]
  l<-length(cf)
  pe<-rmvnorm(1,rep(0,l),cv)
  ncf<-cf+pe
  ncf<-t(ncf)
  nodds<-mdl%*%ncf
  nweights<-exp(nodds)/(1+exp(nodds))
  return(nweights)
}

#Part-2: generating the error to be added to coefficients for the mediator regression
library(mvtnorm)
mcoeffs<-function(mfit){
  coeff<-coef(mfit,matrix=TRUE)
  cv<-vcov(mfit)
  l<-dim(coeff)[1]
  k<-dim(coeff)[2]
  m<-l*k
  #Since the mfit has coefficients corresponding to every level of a multinomial mediator
  #information from the coeff and cv matrix above need to be extracted carefully
  #it is for this reason we at first developed an index, s,
  #this index s now allows to extract the correct set of coefficients and covariances
  #corresponding to a level of a mediator.
  s<-seq(1,m,k)
  #s: 1 9 17 25 33 41 49 57 65 73 81 89 97
  ncoeff<-matrix(0,k,l)
  cf<-coeff[,1]
  cv1<-as.matrix(cv[s,s])
  me<-rmvnorm(1,rep(0,l),cv[s,s])
  ncoeff[1,]<-coeff[,1]+me
  for(i in 1:(k-1)){
    cf<-coeff[,i+1]
    cv1<-as.matrix(cv[s+i,s+i])
    men<-rmvnorm(1,rep(0,l),cv1)
    ncoeff[i+1,]<-cf+men
  }
  ncoeff<-t(ncoeff)
  return(ncoeff)
}


#Part-3 program for computing the predicted probabilities in case of multinomial mediator
trial<-function(mdl,ncoef){
  ro<-dim(mdl)[1]
  co<-dim(mdl)[2]
  k<-ncol(ncoef)
  n<-ro/k
  rsm<-matrix(0,n,k)
  fs<-seq(1,ro,k)
  ss<-seq(1,co,k)
  rsm[,1]<-mdl[fs,ss]%*%ncoef[,1]
  for(i in 1:(k-1)){
    pv<-mdl[fs+i,ss+i]%*%ncoef[,(i+1)]
    rsm[,(i+1)]<-pv
  }
  rsm<-exp(rsm)
  rs<-rowSums(rsm)
  s<-matrix(0,n,k)
  for(i in 1:n){
    s[i,]<-rsm[i,]/(1+rs[i])
  }
  rss<-rowSums(s)
  s<-as.data.frame(s)
  s$V9<-1-rss
  s<-as.matrix(s)
  return(s)
}


# Code for computing the variance
# Once the weights are recomputed from each simulation using the above code then the outcome regression
# needs to be performed for each simulation. Outcome regression coeffcients and the variance-covariance
# matrix from each of this simulation are stored to compute the final variance estimate that accounts for
# varying weights. To do the variance computation using the stored estimates we used the following code.
#Code for exposure regression before conducting simulations
Aregdr<-glm(Atemp~factor(L01)+factor(L02)+factor(L03)
            +factor(L04)+factor(L05),data=mydata,family=binomial)
#code for mediator regression from each time point before simulations
m1reg<-vglm(M1~Atemp+factor(L01)+factor(L02)
            +factor(L03)+factor(L04)+factor(L05)
            ,data=mydata,family=multinomial)
m2reg<-vglm(M2~factor(M1)+Atemp+Y1+factor(L01)
            +factor(L02)+factor(L03)+factor(L04)+factor(L05)
            ,data=mydata,family=multinomial(),maxit=500)
m3reg<-vglm(M3~factor(M2)+factor(M1)+Y2+Y1+Atemp+factor(L01)
            +factor(L02)+factor(L03)+factor(L04)
            +factor(L05),data=mydata,family=multinomial())
m4reg<-vglm(M4~factor(M3)+factor(M2)+factor(M1)+Y2+Y1+Y3+Atemp
            +factor(L01)+factor(L02)+factor(L03)
            +factor(L04)+factor(L05),data=mydata,family=multinomial())
mdl1<-model.matrix(m1reg)
mdl2<-model.matrix(m2reg)
mdl3<-model.matrix(m3reg)
mdl4<-model.matrix(m4reg)
mydata$Atemp<-1-mydata$A
m1creg<-vglm(M1~Atemp+factor(C1)+factor(L01)+factor(L02)+factor(L03)
             +factor(L04)+factor(L05),data=mydata,family=multinomial)
m2creg<-vglm(M2~factor(M1)+Atemp+Y1+factor(L01)+factor(L02)+factor(L03)
             +factor(L04)+factor(L05),data=mydata,family=multinomial(),maxit=500)
m3creg<-vglm(M3~factor(M2)+factor(M1)+Y2+Y1+Atemp+factor(L01)+factor(L02)+factor(L03)
             +factor(L04)+factor(L05),data=mydata,family=multinomial())
m4creg<-vglm(M4~factor(M3)+factor(M2)+factor(M1)+Y2+Y1+Y3+Atemp+factor(L01)
             +factor(L02)+factor(L03)+factor(L04)+factor(L05)
             ,data=mydata,family=multinomial())
mdlc1<-model.matrix(m1creg)
mdlc2<-model.matrix(m2creg)
mdlc3<-model.matrix(m3creg)
mdlc4<-model.matrix(m4creg)

M<-500
t<-matrix(c(0,1,2,3),1,4)
#full code of simulations using above code
for(i in 1:M){
  apdr<-awts(Aregdr)
  wa<-ifelse(mydata$A==1,1/apdr,(1-mydata$A)/(1-apdr)) #weight of A.
  wa_95<-ifelse(wa>=quantile(wa,.95),quantile(wa,.95),wa)
  c1<-mcoeffs(m1reg)
  c2<-mcoeffs(m2reg)
  c3<-mcoeffs(m3reg)
  c4<-mcoeffs(m4reg)
  mt1nr<-trial(mdlc1,c1)[cbind(1:nrow(mydata),mydata$M1)]
  mt2nr<-trial(mdlc2,c2)[cbind(1:nrow(mydata),mydata$M2)]
  mt3nr<-trial(mdlc3,c3)[cbind(1:nrow(mydata),mydata$M3)]
  mt4nr<-trial(mdlc4,c4)[cbind(1:nrow(mydata),mydata$M4)]
  mt1dr<-trial(mdl1,c1)[cbind(1:nrow(mydata),mydata$M1)]
  mt2dr<-trial(mdl2,c2)[cbind(1:nrow(mydata),mydata$M2)]
  mt3dr<-trial(mdl3,c3)[cbind(1:nrow(mydata),mydata$M3)]
  mt4dr<-trial(mdl4,c4)[cbind(1:nrow(mydata),mydata$M4)]
  #creation of new weights
  wM1<-mt1nr/mt1dr
  mydata$wM1<-wM1*wa
  wM2<-(mt2nr/mt2dr)*wM1
  mydata$wM2<-wM2*wa
  wM3<-(mt3nr/mt3dr)*wM2
  mydata$wM3<-wM3*wa
  wM4<-(mt4nr/mt4dr)*wM3
  
  mydata$wM4<-wM4*wa
  #truncated weights
  wM1<-mt1nr/mt1dr
  mydata$wM1_95<-wM1*wa_95
  wM2<-(mt2nr/mt2dr)*wM1
  mydata$wM2_95<-wM2*wa_95
  wM3<-(mt3nr/mt3dr)*wM2
  mydata$wM3_95<-wM3*wa_95
  wM4<-(mt4nr/mt4dr)*wM3
  mydata$wM4_95<-wM4*wa_95
  
  
  #creation of final data for analysis
  gdat1<-cbind(A=mydata$A,Y=mydata$Y1,Astar=mydata$A,T=0,
               id=1:nrow(mydata),W=wa,W_95=wa_95)
  gdat2<-cbind(A=mydata$A,Y=mydata$Y1,Astar=1-mydata$A,T=0,
               id=1:nrow(mydata),W=mydata$wM1,W_95=mydata$wM1_95)
  gdat3<-cbind(A=mydata$A,Y=mydata$Y2,Astar=mydata$A,T=1,
               id=1:nrow(mydata),W=wa,W_95=wa_95)
  gdat4<-cbind(A=mydata$A,Y=mydata$Y2,Astar=1-mydata$A,T=1,
               id=1:nrow(mydata),W=mydata$wM2,W_95=mydata$wM2_95)
  gdat5<-cbind(A=mydata$A,Y=mydata$Y3,Astar=mydata$A,T=2,
               id=1:nrow(mydata),W=wa,W_95=wa_95)
  gdat6<-cbind(A=mydata$A,Y=mydata$Y3,Astar=1-mydata$A,T=2,
               id=1:nrow(mydata),W=mydata$wM3,W_95=mydata$wM3_95)
  gdat7<-cbind(A=mydata$A,Y=mydata$Y4,Astar=mydata$A,T=3,
               id=1:nrow(mydata),W=wa,W_95=wa_95)
  gdat8<-cbind(A=mydata$A,Y=mydata$Y4,Astar=1-mydata$A,T=3,
               id=1:nrow(mydata),W=mydata$wM4,W_95=mydata$wM4_95)
  #combining the data
  newmydata<-as.data.frame(rbind(gdat1,gdat2,gdat3,gdat4
                                 ,gdat5,gdat6,gdat7,gdat8))
  newmydata<-newmydata[order(newmydata$id),]
  newmydata<-newmydata[order(newmydata$id),]
  library(geepack)
  fitsepb<-geeglm(Y~A+Astar+T+I(T*A)+I(T*Astar),id=newmydata$id,
                  corstr="independence",family="binomial",
                  weight=W,data=newmydata,scale.fix=T)
  
  
  #extracting coefficients and variance covariance from each simulation
  cf<-summary(fitsepb)$coefficients
  cv<-vcov(fitsepb)
  #DE matrix of direct effects
  #IDE matrix of indirect effects
  # $VC variance covariance terms required
  DE[i,]<-c(cf[2,1],cf[5,1])
  IDE[i,]<-c(cf[3,1],cf[6,1])
  VC[i,]<-c(cv[2,2],cv[3,3],cv[5,5],cv[6,6],
            cv[2,3],cv[2,5],cv[2,6],cv[3,5]
            ,cv[3,6],cv[5,6])
  #truncated 95
  fitsepb_95<-geeglm(Y~A+Astar+T+I(T*A)+I(T*Astar),id=newmydata$id,
                     corstr="independence",family="binomial",
                     weight=W_95,data=newmydata,scale.fix=T)
  #extracting coefficients and variance covariance matrices
  cf<-summary(fitsepb_95)$coefficients
  cv<-vcov(fitsepb_95)
  #DE direct effect coefficients matrix corresponding to truncated data
  #IDE indirect effect matrix
  #Variance covariances corresponding to truncated data model fit.
  DE_95[i,]<-c(cf[2,1],cf[5,1])
  IDE_95[i,]<-c(cf[3,1],cf[6,1])
  VC_95[i,]<-c(cv[2,2],cv[3,3],cv[5,5],cv[6,6]
               ,cv[2,3],cv[2,5],cv[2,6],
               cv[3,5],cv[3,6],cv[5,6])
}

#compute the final estimate of direct, indirect and total effect for each simulation.
#data for these are from the changing weights
#DE: Stores the coefficients required for computing the direct effect estimates
#IDE: Stores the coefficients required for computing the indirect effect estimates
#TE: total Effects
#DE[,1]: Has the values corresponding to alpha_1 from simulations
#DE[,2]: Has the coefficients corresponding to alpha_4
#IDE[,1]: Has the values corresponding to alpha_2 from simulations
#IDE[,2]: Has the values corresponding to alpha_5
DE_final<-DE[,1]+matrix(DE[,2],M,1)%*%t
IDE_final<-IDE[,1]+matrix(IDE[,2],M,1)%*%t
TE_final<-DE_final+IDE_final
#Step-2
#compute the variances of the DE, IDE and TE for all time points
#data for these are from the variance and covariance extract
#from each simulation with new weights
#VC: Variance covariance estimates stored from each simulation
DEV<-VC[,1]+matrix(VC[,3],M,1)%*%(t^2)+matrix(VC[,6],M,1)%*%(2*t)
IDEV<-VC[,2]+matrix(VC[,4],M,1)%*%(t^2)+matrix(VC[,9],M,1)%*%(2*t)
TEV<-VC[,1]+VC[,2]+matrix(VC[,3],M,1)%*%(t^2)+matrix(VC[,6],M,1)%*%(2*t)
+matrix(VC[,4],M,1)%*%(t^2)+matrix(VC[,9],M,1)%*%(2*t)
+matrix(VC[,5],M,1)%*%(2*t)+matrix(VC[,7],M,1)%*%(2*t)
+matrix(VC[,8],M,1)%*%(2*t)+matrix(VC[,10],M,1)%*%(2*t^2)

#Step-3
#Computing the new variance of the direct, indirect, and total effects
#that account for the weight changes
DEvar<-apply(DE_final,2,var)+apply(DEV,2,mean)
IDEvar<-apply(IDE_final,2,var)+apply(IDEV,2,mean)
TEvar<-apply(TE_final,2,var)+apply(TEV,2,mean)
#Step-4 computing the lower and upper bounds of confidence interval
#which accounts for variance due to changing weights.
#In this step I am using the initial estimates computed from the observed data and
#stored in a matrix labelled "final".
#DEL: Direct effect lower; DEU: Direct effect upper; IDEL: Indirect effect lower;
#IDEU: Indirect effect upper
#TEL: Total effect lower; TEU: Total effect upper.
DEL<-final[1,]-sqrt(DEvar)*1.96
DEU<-final[1,]+sqrt(DEvar)*1.96
IDEL<-final[2,]-sqrt(IDEvar)*1.96
IDEU<-final[2,]+sqrt(IDEvar)*1.96
TEL<-final[3,]-sqrt(TEvar)*1.96
TEU<-final[3,]+sqrt(TEvar)*1.96
#Combining all the estimates to be reported
finalest<-rbind(DEL,final[1,],DEU,IDEL,final[2,],IDEU,TEL,final[3,],TEU)
#the above process is then repeated for the truncated estimates.