#---
#  title: "Figures S2"
#author: ""Colin T. Kremer"
#---

library(ggplot2)
library(dplyr)
library(pbkrtest)
library(lme4)

### Add random effect to analysis of growth rates over time

# - evolution expt paper with MAG
dat1<-read.csv("GrowhRates_evolExp.csv")

# QC filters from Maria's script
dat1<-dat1[!is.na(dat1$nitrate) & !is.na(dat1$rep) & 
               !is.na(dat1$model) & dat1$model!="none",]
dat1$model<-as.character(dat1$model)
dat1$days<-dat1$tiempo.en.dias.desde.el.inicio.del.exp
head(dat1)

##### Run for high nitrate:

dat1.L1<-dat1 %>% filter(nitrate=='L1') %>% filter(days <=200)
dat1.L1$replicate<-factor(dat1.L1$rep)

m0<-lmer(r~(1|replicate),data=dat1.L1)
m2<-lmer(r~days+(1|replicate),data=dat1.L1)
summary(m2)
summary(m0)

# replicate's group means are *very* similar - contributes to a random effect term that is very close to zero, causing the singularity error:
dat1.L1 %>% group_by(replicate) %>% summarise(mean(r))


PBmodcomp(m2,m0)

#Parametric bootstrap test; time: 19.49 sec; samples: 1000 extremes: 7;
#Requested samples: 1000 Used samples: 892 Extremes: 7
#large : r ~ days + (1 | replicate)
#small : r ~ (1 | replicate)
#           stat df  p.value   
#  LRT    7.2846  1 0.006955 **
#  PBtest 7.2846    0.008959 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# over the 200-day interval, no longer consistent with results of AIC comparison...:
AICtab(m0,m2)
#dAIC df
#m0 0.0  3 
#m2 8.7  4

# old version:
#dAIC df
#m2 0.0  4 
#m0 9.1  3 

# Corresponding figure:

pd.dat.L1<-data.frame(days=seq(3,198,1))
#pd.dat.L1<-data.frame(days=seq(1,252,1))
pd.dat.L1$r<-predict(m2,newdata=pd.dat.L1,re.form=NA)

pp.L1<-ggplot(dat1.L1,aes(x=days,y=r))+
  geom_point(aes(shape=replicate),colour="#009E73")+
  geom_line(data=pd.dat.L1,colour="#009E73")+
  geom_vline(xintercept=102,colour='gray',size=1, linetype=2)+
  geom_vline(xintercept=114,colour='gray',size=1, linetype=2)+
  geom_vline(xintercept=165,colour='gray',size=1, linetype=3)+
  geom_vline(xintercept=186,colour='gray',size=1, linetype=3)+
  scale_shape_discrete('Replicate')+
  scale_x_continuous(limits=c(0,200),'Time (days)')+
  scale_y_continuous(limits=c(0,2.5),'Growth rate (1/days)')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text=element_text(size=18, family = "Times New Roman"))



### Now for low nitrogen:

dat1.5<-dat1 %>% filter(nitrate=='5') %>% filter(days <=200)
dat1.5$replicate<-factor(dat1.5$rep)
dat1.5$days2<-dat1.5$days/100


dat1.5$replicate<-factor(dat1.5$replicate)
m0.5<-lmer(r~(1|replicate),data=dat1.5)
m1.5<-lmer(r~days2+(1|replicate),data=dat1.5)
m2.5<-lmer(r~days2+I(days2^2)+(1|replicate),data=dat1.5)
summary(m2.5)
summary(m1.5)
summary(m0.5)


# Test linear term
PBmodcomp(m1.5,m0.5)
# Parametric bootstrap test; time: 19.27 sec; samples: 1000 extremes: 0;
# Requested samples: 1000 Used samples: 888 Extremes: 0
# large : r ~ days2 + (1 | replicate)
# small : r ~ (1 | replicate)
# stat df   p.value    
# LRT    35.188  1 2.994e-09 ***
#   PBtest 35.188     0.001125 ** 

# Test quadratic term alone (against linear)
PBmodcomp(m2.5,m1.5)
#Parametric bootstrap test; time: 19.24 sec; samples: 1000 extremes: 0;
#large : r ~ days + I(days^2) + (1 | replicate)
#small : r ~ days + (1 | replicate)
#          stat df   p.value    
#  LRT    15.31  1 9.121e-05 ***
#  PBtest 15.31     0.000999 ***

# Test quadratic against intercept
PBmodcomp(m2.5,m0.5,	nsim=5000)
# Parametric bootstrap test; time: 110.13 sec; samples: 5000 extremes: 0;
# large : r ~ days2 + I(days2^2) + (1 | replicate)
# small : r ~ (1 | replicate)
         # stat df   p.value    
# LRT    50.498  2 1.083e-11 ***
# PBtest 50.498        2e-04 ***
# ---
# Signif. codes:  
# 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



AICtab(m0.5,m1.5,m2.5)
#dAIC df
#m2.5  0.0 5 
#m1.5  9.7 4 
#m0.5 37.9 3 




# Corresponding figure:

range(dat1.5$days2)
pd.dat.5<-data.frame(days2=seq(0.03,1.98,1/100))
pd.dat.5$r<-predict(m2.5,newdata=pd.dat.5,re.form=NA)
pd.dat.5$days<-pd.dat.5$days2*100

pp.5<-ggplot(dat1.5,aes(x=days,y=r))+
  geom_point(aes(shape=replicate),colour="#E69F00")+
  geom_line(data=pd.dat.5,colour="#E69F00")+
  geom_vline(xintercept=102,colour='gray',size=1, linetype=2)+
  geom_vline(xintercept=114,colour='gray',size=1, linetype=2)+
  geom_vline(xintercept=165,colour='gray',size=1, linetype=3)+
  geom_vline(xintercept=186,colour='gray',size=1, linetype=3)+
  scale_shape_discrete('Replicate')+
  scale_x_continuous(limits=c(0,200),'Time (days)')+
  scale_y_continuous(limits=c(0,2.5),'Growth rate (1/days)')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text=element_text(size=18, family = "Times New Roman"))


##### Combine figures:

library(gridExtra)
grid.arrange(pp.L1,pp.5)

g2<-ggplot(dat1,aes(x=days,y=r))+
  geom_point(aes(shape=factor(rep),colour=nitrate),alpha=0.5,size=3)+
  geom_line(data=pd.dat.5,colour="#E69F00")+
  geom_line(data=pd.dat.L1,colour="#009E73")+
  geom_vline(xintercept=102,colour='gray',size=1, linetype=2)+
  geom_vline(xintercept=114,colour='gray',size=1, linetype=2)+
  geom_vline(xintercept=165,colour='gray',size=1, linetype=3)+
  geom_vline(xintercept=186,colour='gray',size=1, linetype=3)+
  scale_shape_discrete('Replicate')+
  scale_x_continuous(limits=c(0,200),'Time (days)')+
  scale_y_continuous(limits=c(0,2.5),expression("Growth rate (d"^"-1"*")"))+
  scale_colour_manual('Nitrate \nlevel',values=c("#E69F00","#009E73"),
                      breaks=c("5","L1"),
                      labels=c("5","L1"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text=element_text(size=18, family = "Times New Roman"))

png(filename = "GR_Time_31C_v2_Maria.png",units = "mm",res = 300, height = 100, width=300) #I've created a pnj with the specified name, with the spcified size and resolution (usually for ptintting is 300)
g2#call the plot-variable to put it in the png
dev.off()


