# In this script we will re-analyse the data behind the
# recent GARD study - key points we will address are:
# lack of an assessment of the no-interaction assumption
# in the pooled analysis and also the imputed MGMT expression value
# of 5.457623

require(survival)
require(prodlim)
require(ggplot2)
require(tidyverse)

load("data_df.rda")
summary(rsi.all)
# Weird negative OS time
dim(rsi.all)#1574 patients
summary(rsi.all$MGMT_Expression)
# 1330 have missing MGMT expression
1330/1574# 84% missing data
hist(rsi.all$MGMT_Expression, xlab="MGMT Expression",
     main="")
abline(v=5.457623,col=2,lwd=2)
# it doesn't matter but its interesting that the imputed value
# used is not the mean/median - how did they actually come up 
# with this value?

# There is a patient with a negative OS time and one with a 0 value!
rsi.all[rsi.all$Time_OS<=0 & is.na(rsi.all$Time_OS)==F,]
rsi.all<-rsi.all[rsi.all$Time_OS>0|is.na(rsi.all$Time_OS)==T ,]

summary(rsi.all)

# Given the authors claim this top-quality evidence - data QC
# clearly wasn't done

# We shall first explore the association of GARD to OS/RFS and
# assess the no interaction assumption
# Include GARD 
alpha=log(rsi.all$RSI)/(-2) - 0.05*2
rsi.all$gard=rsi.all$n*rsi.all$d*(alpha + 0.05*rsi.all$d)

# assess no. events
sum(rsi.all$Event_OS[rsi.all$Received_RT==1 & rsi.all$Site=="endometrial"],na.rm=T)
sum(rsi.all$Event_OS[rsi.all$Received_RT==1 & rsi.all$Site=="breast_TN"],na.rm=T)# 9 events this study will be removed
sum(rsi.all$Event_OS[rsi.all$Received_RT==1 & rsi.all$Site=="glioma"],na.rm=T)
sum(rsi.all$Event_OS[rsi.all$Received_RT==1 & rsi.all$Site=="pancreas"],na.rm=T)
sum(rsi.all$Event_OS[rsi.all$Received_RT==1 & rsi.all$Site=="lung"],na.rm=T)

m0<-(coxph(Surv(Time_OS,Event_OS)~1+strata(Site),
           data=rsi.all[rsi.all$Received_RT==1 & rsi.all$Site!="breast_TN",]))
m1<-(coxph(Surv(Time_OS,Event_OS)~gard+strata(Site),
                  data=rsi.all[rsi.all$Received_RT==1 & rsi.all$Site!="breast_TN",]))
summary(m1)
# there is no signal!
anova(m1,m0)
m2<-(coxph(Surv(Time_OS,Event_OS)~gard*Site+strata(Site),
           data=rsi.all[rsi.all$Received_RT==1 & rsi.all$Site!="breast_TN",]))
summary(m2)
anova(m1,m2)

# now lets move on to the sham GARD story and the idea that GARD is predictive
# of RT effect
# Lets use the same sham doses first - there is a second debate as to whether
# these are appropriate

rsi.all$Received_RT[rsi.all$Received_RT==1]<-"Yes"
rsi.all$Received_RT[rsi.all$Received_RT==0]<-"No"

rsi.all$TD[rsi.all$Site=="endometrial" & rsi.all$Received_RT=="No"]<-54
rsi.all$TD[rsi.all$Site=="glioma" & rsi.all$Received_RT=="No"]<-60
rsi.all$TD[rsi.all$Site=="pancreas" & rsi.all$Received_RT=="No"]<-50
rsi.all$TD[rsi.all$Site=="breast" & rsi.all$Received_RT=="No"]<-50

rsi.all$n[rsi.all$Site=="endometrial" & rsi.all$Received_RT=="No"]<-27
rsi.all$n[rsi.all$Site=="glioma" & rsi.all$Received_RT=="No"]<-30
rsi.all$n[rsi.all$Site=="pancreas" & rsi.all$Received_RT=="No"]<-25
rsi.all$n[rsi.all$Site=="breast" & rsi.all$Received_RT=="No"]<-25

rsi.all$d[rsi.all$Site=="endometrial" & rsi.all$Received_RT=="No"]<-2
rsi.all$d[rsi.all$Site=="glioma" & rsi.all$Received_RT=="No"]<-2
rsi.all$d[rsi.all$Site=="pancreas" & rsi.all$Received_RT=="No"]<-2
rsi.all$d[rsi.all$Site=="breast" & rsi.all$Received_RT=="No"]<-2

# Include GARD 
alpha=log(rsi.all$RSI)/(-2) - 0.05*2
rsi.all$gard=rsi.all$n*rsi.all$d*(alpha + 0.05*rsi.all$d)

# head and neck, breast_TN and lung are irrelevant in this analysis as is the
# breast NKI cohort

rsi.all<-rsi.all[rsi.all$Site!="HN" & rsi.all$Site!="lung" & 
                   rsi.all$Site!="breast_TN" & rsi.all$Source!="NKI_wboost",]


# we now explore the interaction assumptions
km0<-prodlim(Hist(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="endometrial",])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Survival Fraction",xlim=c(0,12),col=c(1,2),
     marktime=T)
mtext("Endometrial")
survdiff(Surv(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="endometrial",])
#p=0.001


km0<-prodlim(Hist(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="glioma",])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Survival Fraction",xlim=c(0,12),col=c(1,2),
     marktime=T)
mtext("Glioma")
survdiff(Surv(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="glioma",])
#p <0.001
km0<-prodlim(Hist(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="pancreas",])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Survival Fraction",xlim=c(0,12),col=c(1,2),
     marktime=T)
mtext("Pancreas")
survdiff(Surv(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="pancreas",])
# p=0.800

km0<-prodlim(Hist(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="endometrial",])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Survival Fraction",xlim=c(0,12),col=c(1,2),
     marktime=T)
mtext("Endometrial")
survdiff(Surv(Time_OS,Event_OS)~Received_RT,data=rsi.all[rsi.all$Site=="endometrial",])
#p=0.001


km0<-prodlim(Hist(Time,Event)~Received_RT,data=rsi.all[rsi.all$Site=="breast" & rsi.all$Source=="Karolinksa",])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Recurrence-Free Fraction",xlim=c(0,12),col=c(1,2),
     marktime=T)
mtext("Breast - Karolinksa")
survdiff(Surv(Time,Event)~Received_RT,data=rsi.all[rsi.all$Site=="breast" & 
                                                     rsi.all$Source=="Karolinksa",])
#p =0.203
km0<-prodlim(Hist(Time,Event)~Received_RT,data=rsi.all[rsi.all$Site=="breast" & rsi.all$Source=="Erasmus",])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Recurrence-Free Fraction",xlim=c(0,12),col=c(1,2),
     marktime=T)
mtext("Breast - Erasmus")
survdiff(Surv(Time,Event)~Received_RT,data=rsi.all[rsi.all$Site=="breast" & rsi.all$Source=="Erasmus",])
#p =0.060

km0<-prodlim(Hist(Time,Event)~Received_RT,data=rsi.all[rsi.all$Site=="endometrial" ,])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Recurrence-Free Fraction",xlim=c(0,12),col=c(1,2),
     marktime=T)
mtext("Endometrial")
survdiff(Surv(Time,Event)~Received_RT,data=rsi.all[rsi.all$Site=="endometrial",])
#p =0.020

# We can see there is an interaction between RT and outcome 

# Explore OS first

m1<-(coxph(Surv(Time_OS,Event_OS)~gard + Received_RT + strata(Site),
           data=rsi.all))
m1b<-(coxph(Surv(Time_OS,Event_OS)~gard*Received_RT + strata(Site),
           data=rsi.all))
anova(m1,m1b)
# lets now test for the no-interaction assumption i.e. is RT
m1c<-(coxph(Surv(Time_OS,Event_OS)~gard*Received_RT*Site + strata(Site),
            data=rsi.all))
anova(m1b,m1c)
# hmmm so there is an interaction with site
m1d<-(coxph(Surv(Time_OS,Event_OS)~Received_RT*Site + strata(Site),
            data=rsi.all))
anova(m1c,m1d)# p = 0.107


# Explore local progression
m2<-(coxph(Surv(Time,Event)~gard + Received_RT + strata(Source),
           data=rsi.all[rsi.all$Site=="endometrial"|
                          c(rsi.all$Site=="breast" & rsi.all$Source=="Erasmus")|
                          c(rsi.all$Site=="breast" & rsi.all$Source=="Karolinksa"),]))
m2b<-(coxph(Surv(Time,Event)~gard*Received_RT + strata(Site),
            data=rsi.all[rsi.all$Site=="endometrial"|
                           c(rsi.all$Site=="breast" & rsi.all$Source=="Erasmus")|
                           c(rsi.all$Site=="breast" & rsi.all$Source=="Karolinksa"),]))
anova(m2,m2b)# p < 0.001
m2c<-(coxph(Surv(Time,Event)~gard*Received_RT*Site + strata(Site),
            data=rsi.all[rsi.all$Site=="endometrial"|
                           c(rsi.all$Site=="breast" & rsi.all$Source=="Erasmus")|
                           c(rsi.all$Site=="breast" & rsi.all$Source=="Karolinksa"),]))
m2d<-(coxph(Surv(Time,Event)~Received_RT*Site + strata(Site),
            data=rsi.all[rsi.all$Site=="endometrial"|
                           c(rsi.all$Site=="breast" & rsi.all$Source=="Erasmus")|
                           c(rsi.all$Site=="breast" & rsi.all$Source=="Karolinksa"),]))
anova(m2c,m2d)# p = 0.379
