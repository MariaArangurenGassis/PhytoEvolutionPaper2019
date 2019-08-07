

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Set up data for plotting & analyses of Fig. S3-S6, Aranguren-Gassis et al. 2019     ####

# Code by CT Kremer and M Aranguren-Gassis, 2019

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# load packages
library(xlsx)
library(lubridate)
library(dplyr)
library(bbmle)
library(ggplot2)
library(gridExtra)
library(devtools)
library(mleTools)
library(emdbook)

# Double exponential model, re-parameterized to depend explicitly on Topt:
decurve<-function(temp,topt,b1,b2,d0,d2){
  res <- b1*exp(b2*temp) - (d0 + ((b1*b2)/d2)*exp((b2-d2)*topt)*exp(d2*temp))
  res
}

# Function to calculate 95% CI's based on standard error
ci.se<-function(x,na.rm=F){
  if(na.rm){
    x<-na.omit(x)
  }
  1.96*(sd(x)/length(x))
}


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Calculate growth rates - 100 Generations ####

# Note: can skip down to '100 Generation analyses' below

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Import fluorescence data for 34 C time series
GR100<-read.xlsx2("100generations_data.xlsx",1,colClasses=c("character","character","character","character","numeric","character","numeric","numeric"))
GR100<-GR100[,1:8]

# focus on 34 C:
GR100<-GR100[GR100$TEMPERATURE==34,]

#converting date strings in real date and then in seconds and then in days
GR100$DATE.AND.TIME<-as.POSIXct(GR100$DATE.AND.TIME,format="%m/%d/%Y %H:%M",tz="EST")
GR100$DATE.AND.TIME<-as.numeric(GR100$DATE.AND.TIME)/86400

#calculating log biomass
GR100$ln.fluor<-log(GR100$Fluorescence)
head(GR100)

# define a function for automating growth rate calculation over first 4 days:
get.gr.rate<-function(xs,ys){
  dat<-data.frame(xs,ys)
  dat<-dat[order(dat$xs),]
  dat<-dat[1:4,]
  
  m1<-lm(ys~xs,data=dat)
  coef(m1)[2]
}

# use get.gr.rate() to calculate growth rate over first 4 observations 
rates100 <- GR100 %>% group_by(EVOLUTION.EXPERIMENT.STRAIN,TEMPERATURE,NITRATE.CONCENTRATION,Wellplate.replicate) %>% summarise(Growth.rate=get.gr.rate(DATE.AND.TIME,ln.fluor))

# format the result
# library(plyr)
rates100<-plyr::rename(rates100,c('TEMPERATURE'='Temperature','Wellplate.replicate'='Well.Replicate','EVOLUTION.EXPERIMENT.STRAIN'='Evol.Strain'))
rates100<-rates100[,-3]
head(rates100)

# Save the growth rates at 34ºC
# write.xlsx2(rates100, file="100g_TPCs_GrowthRates_34_firstDays.xlsx",sheetName = "Sheet1",
#            col.names = TRUE, row.names = TRUE, append = FALSE) 

#### Load previously calculated growth rates:

# specific to 34ºC
rates100<-read.xlsx2(file="100g_TPCs_GrowthRates_34_firstDays.xlsx",sheetName = "Sheet1",colClasses=c("character","character","numeric","character","numeric"))
rates100<-rates100[,-1]

# for all temperatures (drop 34 C later and combine with rates100)
tmp100<-read.xlsx2(file="Curve100G_Chsim_NoPre.xlsx",1,colClasses=c("character","numeric","numeric","numeric","numeric"))
tmp100<-tmp100[,1:3]

# Decompose information contained in ID into separate columns
tmp100$ID<-as.character(tmp100$ID)
a<-strsplit(tmp100$ID," +")
tmp100$Evol.Strain<-unlist(lapply(a, '[[', 1))
tmp100$Temperature<-unlist(lapply(a, '[[', 2))
tmp100$Well.Replicate<-unlist(lapply(a, '[[', 4))

# drop values from 34 C
tmp100$Temperature<-as.numeric(tmp100$Temperature)
tmp100<-tmp100 %>% filter(Temperature!=34)
head(tmp100)

### Combine these two different sets of rates:

# thin columns of tmp100
tmp100<-tmp100[,names(rates100)]

# take 34 C values from rates100 and combine with tmp100
tmp100<-rbind(tmp100,rates100)
head(tmp100)

# save output:
#write.csv(tmp100,'growth_rates_at_100gen_042519.csv',row.names=F)


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### 100 Generation analysis ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Load growth rate data:
tmp100<-read.csv('growth_rates_at_100gen_042519.csv')
head(tmp100)

# split data sets for separate analyses
tmp100.L1<-tmp100 %>% filter(Evol.Strain %in% c('Collection','L1-1','L1-2','L1-3','L1-4'))
tmp100.5<-tmp100 %>% filter(Evol.Strain %in% c('Collection','5-1','5-2','5-3','5-4'))


## Construct dummy variable columns for model fitting

# for evolved vs. control comparison
tmp100.5$evolved.Q<-ifelse(grepl(pattern = "5",x = tmp100.5$Evol.Strain),1,0)
tmp100.L1$evolved.Q<-ifelse(grepl(pattern = "L1",x = tmp100.L1$Evol.Strain),1,0)

# for single population comparisons
tmp100.5$l1.Q<-ifelse(grepl(pattern = "5-1",x = tmp100.5$Evol.Strain),1,0)
tmp100.5$l2.Q<-ifelse(grepl(pattern = "5-2",x = tmp100.5$Evol.Strain),1,0)
tmp100.5$l3.Q<-ifelse(grepl(pattern = "5-3",x = tmp100.5$Evol.Strain),1,0)
tmp100.5$l4.Q<-ifelse(grepl(pattern = "5-4",x = tmp100.5$Evol.Strain),1,0)
tmp100.L1$l1.Q<-ifelse(grepl(pattern = "L1-1",x = tmp100.L1$Evol.Strain),1,0)
tmp100.L1$l2.Q<-ifelse(grepl(pattern = "L1-2",x = tmp100.L1$Evol.Strain),1,0)
tmp100.L1$l3.Q<-ifelse(grepl(pattern = "L1-3",x = tmp100.L1$Evol.Strain),1,0)
tmp100.L1$l4.Q<-ifelse(grepl(pattern = "L1-4",x = tmp100.L1$Evol.Strain),1,0)


#### (1) N-limited results ####

# allow each population to have a unique TPC

topt.guess<-28

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,topt.l1=0,topt.l2=0,topt.l3=0,topt.l4=0,
            b1=NA,b1.l1=0,b1.l2=0,b1.l3=0,b1.l4=0,
            b2=NA,b2.l1=0,b2.l2=0,b2.l3=0,b2.l4=0,
            d0=NA,d0.l1=0,d0.l2=0,d0.l3=0,d0.l4=0,
            d2=NA,d2.l1=0,d2.l2=0,d2.l3=0,d2.l4=0,
            s=log(2))

# run fit:
fit2<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,
                                                         topt+topt.l1*l1.Q+topt.l2*l2.Q+topt.l3*l3.Q+topt.l4*l4.Q,
                                                         exp(b1+b1.l1*l1.Q+b1.l2*l2.Q+b1.l3*l3.Q+b1.l4*l4.Q),
                                                         exp(b2+b2.l1*l1.Q+b2.l2*l2.Q+b2.l3*l3.Q+b2.l4*l4.Q),
                                                         exp(d0+d0.l1*l1.Q+d0.l2*l2.Q+d0.l3*l3.Q+d0.l4*l4.Q),
                                                         exp(d2+d2.l1*l1.Q+d2.l2*l2.Q+d2.l3*l3.Q+d2.l4*l4.Q)),
                                            sd=exp(s)),grids=grids,start=start,data=tmp100.5,
                control=list(maxit=10000))
head(tmp100.5)

# extract parameters for polished fit
cfg2<-as.list(coef(fit2$res.best))
guesses2<-list(topt=cfg2$topt,topt.l1=cfg2$topt.l1,topt.l2=cfg2$topt.l2,topt.l3=cfg2$topt.l3,topt.l4=cfg2$topt.l4,
               b1=exp(cfg2$b1),b1.l1=exp(cfg2$b1+cfg2$b1.l1)-exp(cfg2$b1),b1.l2=exp(cfg2$b1+cfg2$b1.l2)-exp(cfg2$b1),b1.l3=exp(cfg2$b1+cfg2$b1.l3)-exp(cfg2$b1),b1.l4=exp(cfg2$b1+cfg2$b1.l4)-exp(cfg2$b1),
               b2=exp(cfg2$b2),b2.l1=exp(cfg2$b2+cfg2$b2.l1)-exp(cfg2$b2),b2.l2=exp(cfg2$b2+cfg2$b2.l2)-exp(cfg2$b2),b2.l3=exp(cfg2$b2+cfg2$b2.l3)-exp(cfg2$b2),b2.l4=exp(cfg2$b2+cfg2$b2.l4)-exp(cfg2$b2),
               d0=exp(cfg2$d0),d0.l1=exp(cfg2$d0+cfg2$d0.l1)-exp(cfg2$d0),d0.l2=exp(cfg2$d0+cfg2$d0.l2)-exp(cfg2$d0),d0.l3=exp(cfg2$d0+cfg2$d0.l3)-exp(cfg2$d0),d0.l4=exp(cfg2$d0+cfg2$d0.l4)-exp(cfg2$d0),
               d2=exp(cfg2$d2),d2.l1=exp(cfg2$d2+cfg2$d2.l1)-exp(cfg2$d2),d2.l2=exp(cfg2$d2+cfg2$d2.l2)-exp(cfg2$d2),d2.l3=exp(cfg2$d2+cfg2$d2.l3)-exp(cfg2$d2),d2.l4=exp(cfg2$d2+cfg2$d2.l4)-exp(cfg2$d2),
               s=exp(cfg2$s))

# polish best fit model using formula interface
fit2P.5<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.l1*l1.Q+topt.l2*l2.Q+topt.l3*l3.Q+topt.l4*l4.Q,
                                             b1+b1.l1*l1.Q+b1.l2*l2.Q+b1.l3*l3.Q+b1.l4*l4.Q,
                                             b2+b2.l1*l1.Q+b2.l2*l2.Q+b2.l3*l3.Q+b2.l4*l4.Q,
                                             d0+d0.l1*l1.Q+d0.l2*l2.Q+d0.l3*l3.Q+d0.l4*l4.Q,
                                             d2+d2.l1*l1.Q+d2.l2*l2.Q+d2.l3*l3.Q+d2.l4*l4.Q),
                                sd=s),
              start=guesses2,data=tmp100.5,control=list(maxit=0))
summary(fit2P.5)

# extract predicted curves for plotting:
pd.df2.5<-merge(expand.grid(Evol.Strain=unique(tmp100.5$Evol.Strain),Temperature=seq(10,34,0.1)),unique(tmp100.5[,c("Evol.Strain","evolved.Q","l1.Q","l2.Q","l3.Q","l4.Q")]))
pd.df2.5<-pd.df2.5[order(pd.df2.5$Evol.Strain,pd.df2.5$Temperature),]
pd.df2.5$Growth.rate<-predict(fit2P.5,newdata=pd.df2.5)
head(pd.df2.5)

# aggregate observations to avoid noisy plot
tmp100.5.avg<-tmp100.5 %>% group_by(Evol.Strain,Temperature,evolved.Q) %>% summarise(mn.Growth.rate=mean(Growth.rate),se.Growth.rate=ci.se(Growth.rate))

#### .......Create Fig. 1, panel a ####
F1.a<-ggplot(tmp100.5,aes(y=Growth.rate))+
  geom_point(data=tmp100.5.avg,aes(x=Temperature,y=mn.Growth.rate,
                                    colour=factor(evolved.Q),shape=factor(evolved.Q)))+
  geom_errorbar(data=tmp100.5.avg,aes(x=Temperature,y=mn.Growth.rate,
                                       ymin=mn.Growth.rate-se.Growth.rate,
                                       ymax=mn.Growth.rate+se.Growth.rate,
                                       colour=factor(evolved.Q)),width=0.3,
                show.legend=F)+
  geom_line(data=pd.df2.5,aes(x=Temperature,colour=factor(evolved.Q),linetype=factor(evolved.Q),group=Evol.Strain))+
  geom_hline(yintercept = 0)+
  scale_shape_manual('Treatments',values=c(4,1),
                     breaks=c("0","1"),labels = c("Control","N-limited"))+  
  scale_linetype_manual('Treatments',values=c(2,1),
                        breaks=c("0","1"),labels = c("Control","N-limited"))+
  scale_colour_manual('Treatments',values=c("black","#E69F00"),
                      breaks=c("0","1"),labels = c("Control","N-limited"))+
  scale_x_continuous('Temperature (C)',limits=c(10,35))+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('a) N-Limited')


#### (2) N-replete results ####

# allow each population to have a unique TPC

topt.guess<-28

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,topt.l1=0,topt.l2=0,topt.l3=0,topt.l4=0,
            b1=NA,b1.l1=0,b1.l2=0,b1.l3=0,b1.l4=0,
            b2=NA,b2.l1=0,b2.l2=0,b2.l3=0,b2.l4=0,
            d0=NA,d0.l1=0,d0.l2=0,d0.l3=0,d0.l4=0,
            d2=NA,d2.l1=0,d2.l2=0,d2.l3=0,d2.l4=0,
            s=log(2))

# run fit:
fit2<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,
                                                         topt+topt.l1*l1.Q+topt.l2*l2.Q+topt.l3*l3.Q+topt.l4*l4.Q,
                                                         exp(b1+b1.l1*l1.Q+b1.l2*l2.Q+b1.l3*l3.Q+b1.l4*l4.Q),
                                                         exp(b2+b2.l1*l1.Q+b2.l2*l2.Q+b2.l3*l3.Q+b2.l4*l4.Q),
                                                         exp(d0+d0.l1*l1.Q+d0.l2*l2.Q+d0.l3*l3.Q+d0.l4*l4.Q),
                                                         exp(d2+d2.l1*l1.Q+d2.l2*l2.Q+d2.l3*l3.Q+d2.l4*l4.Q)),
                                            sd=exp(s)),grids=grids,start=start,data=tmp100.L1,
                control=list(maxit=10000))

# extract parameters for polished fit
cfg2<-as.list(coef(fit2$res.best))
guesses2<-list(topt=cfg2$topt,topt.l1=cfg2$topt.l1,topt.l2=cfg2$topt.l2,topt.l3=cfg2$topt.l3,topt.l4=cfg2$topt.l4,
               b1=exp(cfg2$b1),b1.l1=exp(cfg2$b1+cfg2$b1.l1)-exp(cfg2$b1),b1.l2=exp(cfg2$b1+cfg2$b1.l2)-exp(cfg2$b1),b1.l3=exp(cfg2$b1+cfg2$b1.l3)-exp(cfg2$b1),b1.l4=exp(cfg2$b1+cfg2$b1.l4)-exp(cfg2$b1),
               b2=exp(cfg2$b2),b2.l1=exp(cfg2$b2+cfg2$b2.l1)-exp(cfg2$b2),b2.l2=exp(cfg2$b2+cfg2$b2.l2)-exp(cfg2$b2),b2.l3=exp(cfg2$b2+cfg2$b2.l3)-exp(cfg2$b2),b2.l4=exp(cfg2$b2+cfg2$b2.l4)-exp(cfg2$b2),
               d0=exp(cfg2$d0),d0.l1=exp(cfg2$d0+cfg2$d0.l1)-exp(cfg2$d0),d0.l2=exp(cfg2$d0+cfg2$d0.l2)-exp(cfg2$d0),d0.l3=exp(cfg2$d0+cfg2$d0.l3)-exp(cfg2$d0),d0.l4=exp(cfg2$d0+cfg2$d0.l4)-exp(cfg2$d0),
               d2=exp(cfg2$d2),d2.l1=exp(cfg2$d2+cfg2$d2.l1)-exp(cfg2$d2),d2.l2=exp(cfg2$d2+cfg2$d2.l2)-exp(cfg2$d2),d2.l3=exp(cfg2$d2+cfg2$d2.l3)-exp(cfg2$d2),d2.l4=exp(cfg2$d2+cfg2$d2.l4)-exp(cfg2$d2),
               s=exp(cfg2$s))

# polish best fit model using formula interface
fit2P.L1<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.l1*l1.Q+topt.l2*l2.Q+topt.l3*l3.Q+topt.l4*l4.Q,
                                              b1+b1.l1*l1.Q+b1.l2*l2.Q+b1.l3*l3.Q+b1.l4*l4.Q,
                                              b2+b2.l1*l1.Q+b2.l2*l2.Q+b2.l3*l3.Q+b2.l4*l4.Q,
                                              d0+d0.l1*l1.Q+d0.l2*l2.Q+d0.l3*l3.Q+d0.l4*l4.Q,
                                              d2+d2.l1*l1.Q+d2.l2*l2.Q+d2.l3*l3.Q+d2.l4*l4.Q),
                                 sd=s),
               start=guesses2,data=tmp100.L1,control=list(maxit=0))
summary(fit2P.L1)

# extract predicted curves for plotting:
pd.df2.L1<-merge(expand.grid(Evol.Strain=unique(tmp100.L1$Evol.Strain),Temperature=seq(10,35,0.1)),unique(tmp100.L1[,c("Evol.Strain","evolved.Q","l1.Q","l2.Q","l3.Q","l4.Q")]))
pd.df2.L1<-pd.df2.L1[order(pd.df2.L1$Evol.Strain,pd.df2.L1$Temperature),]
pd.df2.L1$Growth.rate<-predict(fit2P.L1,newdata=pd.df2.L1)

# aggregate observations to avoid noisy plot
tmp100.L1.avg<-tmp100.L1 %>% group_by(Evol.Strain,Temperature,evolved.Q) %>% summarise(mn.Growth.rate=mean(Growth.rate),se.Growth.rate=ci.se(Growth.rate))

#### .......Create Fig. 1, panel b ####
F1.b<-ggplot(tmp100.L1,aes(y=Growth.rate))+
  geom_point(data=tmp100.L1.avg,aes(x=Temperature,y=mn.Growth.rate,
                                    colour=factor(evolved.Q),shape=factor(evolved.Q)))+
  geom_errorbar(data=tmp100.L1.avg,aes(x=Temperature,y=mn.Growth.rate,
                                       ymin=mn.Growth.rate-se.Growth.rate,
                                       ymax=mn.Growth.rate+se.Growth.rate,
                                       colour=factor(evolved.Q)),width=0.3,
                show.legend = F)+
  geom_line(data=pd.df2.L1,aes(x=Temperature,colour=factor(evolved.Q),linetype=factor(evolved.Q),group=Evol.Strain))+
  geom_hline(yintercept = 0)+
  scale_shape_manual('Treatments',values=c(4,1),
                     breaks=c("0","1"),labels = c("Control","N-Replete"))+  
  scale_linetype_manual('Treatments',values=c(2,1),
                        breaks=c("0","1"),labels = c("Control","N-Replete"))+
  scale_colour_manual('Treatments',values=c("black","#009E73"),
                      breaks=c("0","1"),labels = c("Control","N-Replete"))+
  scale_x_continuous('Temperature (C)',limits=c(10,35))+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('b) N-Replete')


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Fig. S3: Pairwise TPC comparison (high N evolved vs. control)  ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Approach:
# - pull out data for each individual evolved line,
# - pool it with data on the control line
# - fit (a) a single TPC and (b) 2 TPCs to this data, 
#   either ignoring or acknowledging differences between
#   evolved and control lines
# - test which model fits better
# - plot results

# Here are the IDs of the high N evolved lines
ids.L1<-sort(unique(tmp100.L1$Evol.Strain[tmp100.L1$Evol.Strain!='Collection']))

topt.guess<-28
f0.L1.aicc<-c()
f1.L1.aicc<-c()
f0s.L1<-list()
f1s.L1<-list()
plot.L1.list<-list()
for(i in 1:length(ids.L1)){
  tmpdat<-tmp100.L1[tmp100.L1$Evol.Strain %in% c('Collection',as.character(ids.L1[i])),]
  
  ## Fit single curve
  
  # set up a search grid of parameter guesses
  grids<-list(b1=log(seq(0.01,5,0.5)),b2=log(seq(0.1,0.5,0.2)),
              d0=log(seq(0.01,5,0.5)),d2=log(seq(0.1,0.7,0.2)))
  start<-list(topt=topt.guess,b1=NA,b2=NA,d0=NA,d2=NA,s=log(5))
  
  # run fit:
  fit0.L1<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt,exp(b1),exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tmpdat,control=list(maxit=10000))
  
  # extract parameters for polished fit
  cfg<-as.list(coef(fit0.L1$res.best))
  guesses<-list(topt=cfg$topt,b1=exp(cfg$b1),b2=exp(cfg$b2),d0=exp(cfg$d0),d2=exp(cfg$d2),s=exp(cfg$s))
  
  # polish best fit model using formula interface
  fit0.L1.P<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
                 start=guesses,data=tmpdat,control=list(maxit=0))
  summary(fit0.L1.P)
  
  ## Fit 2-curve
  
  # set up search of a grid of parameter guesses
  grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
              d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
  start<-list(topt=topt.guess,topt.T=0,b1=NA,b1.T=0,b2=NA,b2.T=0,d0=NA,d0.T=0,d2=NA,d2.T=0,s=log(2))
  
  # run fit:
  fit1.L1<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,exp(b1+b1.T*evolved.Q),exp(b2+b2.T*evolved.Q),exp(d0+d0.T*evolved.Q),exp(d2+d2.T*evolved.Q)),sd=exp(s)),grids=grids,start=start,data=tmpdat,control=list(maxit=10000))
  
  # extract parameters for polished fit
  cfg1<-as.list(coef(fit1.L1$res.best))
  guesses1<-list(topt=cfg1$topt,topt.T=cfg1$topt.T,
                 b1=exp(cfg1$b1),b1.T=exp(cfg1$b1+cfg1$b1.T)-exp(cfg1$b1),
                 b2=exp(cfg1$b2),b2.T=exp(cfg1$b2+cfg1$b2.T)-exp(cfg1$b2),
                 d0=exp(cfg1$d0),d0.T=exp(cfg1$d0+cfg1$d0.T)-exp(cfg1$d0),
                 d2=exp(cfg1$d2),d2.T=exp(cfg1$d2+cfg1$d2.T)-exp(cfg1$d2),
                 s=exp(cfg1$s))
  
  # polish best fit model using formula interface
  fit1.L1.P<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,b1+b1.T*evolved.Q,b2+b2.T*evolved.Q,d0+d0.T*evolved.Q,d2+d2.T*evolved.Q),sd=s),
                 start=guesses1,data=tmpdat,control=list(maxit=0))
  
  # Compare results:
  print(as.character(ids.L1[i]))
  print(AICctab(fit0.L1.P,fit1.L1.P,nobs=nrow(tmpdat)))
  
  # Save results:
  f0s.L1[[i]]<-fit0.L1.P
  f1s.L1[[i]]<-fit1.L1.P
  f0.L1.aicc<-append(f0.L1.aicc,AICc(fit0.L1.P))
  f1.L1.aicc<-append(f1.L1.aicc,AICc(fit1.L1.P))
  
  # Generate figure:
  model<-f1s.L1[[i]]

  # predict curves given model fit
  pd.df1.L1<-merge(expand.grid(Evol.Strain=unique(tmpdat$Evol.Strain),Temperature=seq(10,34,0.1)),unique(tmpdat[,c("Evol.Strain","evolved.Q")]))
  pd.df1.L1<-pd.df1.L1[order(pd.df1.L1$Evol.Strain,pd.df1.L1$Temperature),]
  pd.df1.L1$Growth.rate<-predict(model,newdata=pd.df1.L1)
  head(pd.df1.L1)
  
  # Repeat fit of Control population, to obtain mle2 object based only on Control data:
  # but force maxit=0 to prevent actual change to parameters.
  cfs<-coef(f1s.L1[[i]])[c('topt','b1','b2','d0','d2','s')]
  fitC<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
             start=as.list(cfs),data=tmp100.L1[tmp100.L1$Evol.Strain=='Collection',],control=list(maxit=0))
  
  # Likewise for evolved population:
  cfs2<-coef(f1s.L1[[i]])[c('topt.T','b1.T','b2.T','d0.T','d2.T','s')]
  cfsE<-cfs+cfs2
  cfsE[names(cfsE)=='s']<-coef(f1s.L1[[i]])['s'][[1]] # note: allowing this to float from comparison to comparison causes slight changes in the apparent confidence band around the control population.
  fitE<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
             start=as.list(cfsE),data=tmp100.L1[tmp100.L1$Evol.Strain==as.character(ids.L1[i]),],control=list(maxit=0))
  
  # check that parameters haven't changed...
  coef(fitE)
  cfsE
  
  # Now, generate confidence bands making use of these updated/separated fits:
  xs<-seq(10,34,0.1)
  cfsC<-coef(fitC)
  cfsE<-coef(fitE)
  dvsC<-deltavar(fun=decurve(xs,topt,b1,b2,d0,d2),meanval=cfsC,Sigma=vcov(fitC))
  dvsE<-deltavar(fun=decurve(xs,topt,b1,b2,d0,d2),meanval=cfsE,Sigma=vcov(fitE))
  bds<-data.frame(Temperature=c(xs,xs),se=c(sqrt(dvsC),sqrt(dvsE)),Evol.Strain=c(rep('Collection',length(xs)),rep(as.character(ids.L1[i]),length(xs))))
  pds<-merge(pd.df1.L1,bds)
  head(pds)
  
  # visualize fits:
  plot.L1<-ggplot(tmpdat,aes(x=Temperature,y=Growth.rate))+
    geom_ribbon(data=pds,aes(ymin=Growth.rate-1.96*se,ymax=Growth.rate+1.96*se,group=Evol.Strain,fill=Evol.Strain),alpha=0.2)+
    geom_point(aes(colour=factor(evolved.Q),shape=Evol.Strain))+
    geom_line(data=pd.df1.L1,aes(colour=factor(evolved.Q)))+
    geom_hline(yintercept = 0)+
    scale_shape_discrete('Lineage')+
    scale_linetype_discrete('Lineage')+
    scale_colour_manual('Lineage',values=c("black",'#009E73'),
                        breaks=c("1","0"),labels=c(as.character(ids.L1[i]),'Collection'))+
    scale_fill_manual('Lineage',values=c(gray(0.2),'#009E73'))+
    scale_x_continuous('Temperature (C)')+
    scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
    theme_bw()+
    ggtitle(gsub(pattern = 'L1-',x = ids.L1[i],replacement = 'Population '))
  print(plot.L1)
  plot.L1.list[[i]]<-plot.L1
  
  ggsave(paste("plot_",ids.L1[i],"_test.pdf",sep=''),plot.L1)
}

# 7/10/19 - updated results:

# [1] "L1-1"
# dAICc df
# fit1.L1.P  0.0  11
# fit0.L1.P 18.2  6 
# 
# [1] "L1-2"
# dAICc df
# fit0.L1.P  0.0  6 
# fit1.L1.P  0.8  11
# 
# [1] "L1-3"
# dAICc df
# fit0.L1.P  0.0  6 
# fit1.L1.P  7.9  11
# 
# [1] "L1-4"
# dAICc df
# fit1.L1.P  0.0  11
# fit0.L1.P 24.3  6 

# remove legend
plot.L1.list[[1]]<-plot.L1.list[[1]]+theme(legend.position = 'none')
plot.L1.list[[2]]<-plot.L1.list[[2]]+theme(legend.position = 'none')
plot.L1.list[[3]]<-plot.L1.list[[3]]+theme(legend.position = 'none')
plot.L1.list[[4]]<-plot.L1.list[[4]]+theme(legend.position = 'none')

# export multi-panel figure
grid.arrange(plot.L1.list[[1]],plot.L1.list[[2]],plot.L1.list[[3]],plot.L1.list[[4]],nrow=2)
FS3<-arrangeGrob(plot.L1.list[[1]],plot.L1.list[[2]],plot.L1.list[[3]],plot.L1.list[[4]],nrow=2)

ggsave(paste("Figure_S3.pdf",sep=''),FS3,width=8.25,height=7.75)



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Fig. S4: Pairwise TPC comparison (low N evolved vs. control)  ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Approach:
# - pull out data for each individual evolved line,
# - pool it with data on the control line
# - fit (a) a single TPC and (b) 2 TPCs to this data, 
#   either ignoring or acknowledging differences between
#   evolved and control lines
# - test which model fits better
# - plot results

# Here are the IDs of the low N evolved lines
ids<-sort(unique(tmp100.5$Evol.Strain[tmp100.5$Evol.Strain!='Collection']))

topt.guess<-28
f0.aicc<-c()
f1.aicc<-c()
f0s<-list()
f1s<-list()
plot.5.list<-list()
for(i in 1:length(ids)){
  tmpdat<-tmp100.5[tmp100.5$Evol.Strain %in% c('Collection',as.character(ids[i])),]
  
  ## Fit single curve
  
  # set up a search grid of parameter guesses
  grids<-list(b1=log(seq(0.01,5,0.5)),b2=log(seq(0.1,0.5,0.2)),
              d0=log(seq(0.01,5,0.5)),d2=log(seq(0.1,0.7,0.2)))
  start<-list(topt=topt.guess,b1=NA,b2=NA,d0=NA,d2=NA,s=log(5))
  
  # run fit:
  fit0.5<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt,exp(b1),exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tmpdat,control=list(maxit=10000))
  
  # extract parameters for polished fit
  cfg<-as.list(coef(fit0.5$res.best))
  guesses<-list(topt=cfg$topt,b1=exp(cfg$b1),b2=exp(cfg$b2),d0=exp(cfg$d0),d2=exp(cfg$d2),s=exp(cfg$s))
  
  # polish best fit model using formula interface
  fit0.5.P<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
              start=guesses,data=tmpdat,control=list(maxit=0))
  summary(fit0.5.P)
  
  ## Fit 2-curve
  
  # set up search of a grid of parameter guesses
  grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
              d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
  start<-list(topt=topt.guess,topt.T=0,b1=NA,b1.T=0,b2=NA,b2.T=0,d0=NA,d0.T=0,d2=NA,d2.T=0,s=log(2))
  
  # run fit:
  fit1.5<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,exp(b1+b1.T*evolved.Q),exp(b2+b2.T*evolved.Q),exp(d0+d0.T*evolved.Q),exp(d2+d2.T*evolved.Q)),sd=exp(s)),grids=grids,start=start,data=tmpdat,control=list(maxit=10000))
  
  # extract parameters for polished fit
  cfg1<-as.list(coef(fit1.5$res.best))
  guesses1<-list(topt=cfg1$topt,topt.T=cfg1$topt.T,
                 b1=exp(cfg1$b1),b1.T=exp(cfg1$b1+cfg1$b1.T)-exp(cfg1$b1),
                 b2=exp(cfg1$b2),b2.T=exp(cfg1$b2+cfg1$b2.T)-exp(cfg1$b2),
                 d0=exp(cfg1$d0),d0.T=exp(cfg1$d0+cfg1$d0.T)-exp(cfg1$d0),
                 d2=exp(cfg1$d2),d2.T=exp(cfg1$d2+cfg1$d2.T)-exp(cfg1$d2),
                 s=exp(cfg1$s))
  
  # polish best fit model using formula interface
  fit1.5.P<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,b1+b1.T*evolved.Q,b2+b2.T*evolved.Q,d0+d0.T*evolved.Q,d2+d2.T*evolved.Q),sd=s),
                start=guesses1,data=tmpdat,control=list(maxit=0))
  summary(fit1.5.P)
  
  
  # Compare results:
  print(as.character(ids[i]))
  print(AICctab(fit0.5.P,fit1.5.P,nobs=nrow(tmpdat)))
  
  # Save results:
  f0s[[i]]<-fit0.5.P
  f1s[[i]]<-fit1.5.P
  f0.aicc<-append(f0.aicc,AICc(fit0.5.P))
  f1.aicc<-append(f1.aicc,AICc(fit1.5.P))
  
  # Generate figure:
  model<-f1s[[i]]
  #tmpdat<-tmp100.5[tmp100.5$Evol.Strain %in% c('Collection',as.character(ids[i])),]
  
  # predict curves given model fit
  pd.df1.5<-merge(expand.grid(Evol.Strain=unique(tmpdat$Evol.Strain),Temperature=seq(10,34,0.1)),unique(tmpdat[,c("Evol.Strain","evolved.Q")]))
  pd.df1.5<-pd.df1.5[order(pd.df1.5$Evol.Strain,pd.df1.5$Temperature),]
  pd.df1.5$Growth.rate<-predict(model,newdata=pd.df1.5)
  head(pd.df1.5)
  
  # Repeat fit of Control population, to obtain mle2 object based only on Control data:
  # but force maxit=0 to prevent actual change to parameters.
  cfs<-coef(f1s[[i]])[c('topt','b1','b2','d0','d2','s')]
  fitC<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
             start=as.list(cfs),data=tmp100.5[tmp100.5$Evol.Strain=='Collection',],control=list(maxit=0))
  
  # Likewise for evolved population:
  cfs2<-coef(f1s[[i]])[c('topt.T','b1.T','b2.T','d0.T','d2.T','s')]
  cfsE<-cfs+cfs2
  cfsE[names(cfsE)=='s']<-coef(f1s[[i]])['s'][[1]] # note: allowing this to float from comparison to comparison causes slight changes in the apparent confidence band around the control population.
  fitE<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
             start=as.list(cfsE),data=tmp100.5[tmp100.5$Evol.Strain==as.character(ids[i]),],control=list(maxit=0))
  
  # check that parameters haven't changed...
  coef(fitE)
  cfsE
  
  # Now, generate confidence bands making use of these updated/separated fits:
  xs<-seq(10,34,0.1)
  cfsC<-coef(fitC)
  cfsE<-coef(fitE)
  dvsC<-deltavar(fun=decurve(xs,topt,b1,b2,d0,d2),meanval=cfsC,Sigma=vcov(fitC))
  dvsE<-deltavar(fun=decurve(xs,topt,b1,b2,d0,d2),meanval=cfsE,Sigma=vcov(fitE))
  bds<-data.frame(Temperature=c(xs,xs),se=c(sqrt(dvsC),sqrt(dvsE)),Evol.Strain=c(rep('Collection',length(xs)),rep(as.character(ids[i]),length(xs))))
  pds<-merge(pd.df1.5,bds)
  head(pds)
  
  # visualize fits:
  plot.5<-ggplot(tmpdat,aes(x=Temperature,y=Growth.rate))+
    geom_ribbon(data=pds,aes(ymin=Growth.rate-1.96*se,ymax=Growth.rate+1.96*se,group=Evol.Strain,fill=Evol.Strain),alpha=0.2)+
    geom_point(aes(colour=factor(evolved.Q),shape=Evol.Strain))+
    geom_line(data=pd.df1.5,aes(colour=factor(evolved.Q)))+
    geom_hline(yintercept = 0)+
    scale_shape_discrete('Lineage')+
    scale_linetype_discrete('Lineage')+
    scale_colour_manual('Lineage',values=c("black","#E69F00"),
                        breaks=c("1","0"),labels=c(as.character(ids[i]),'Collection'))+
    scale_fill_manual('Lineage',values=c("#E69F00",gray(0.2)))+
    scale_x_continuous('Temperature (C)')+
    scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
    theme_bw()+
    ggtitle(gsub(pattern = '5-',x = ids[i],replacement = 'Population '))
  print(plot.5)
  plot.5.list[[i]]<-plot.5
  
  ggsave(paste("/Users/colin/Research/Active/plot_",ids[i],"_test.pdf",sep=''),plot.5)
}


# 7/10/19
# in 3 out of 4 populations, the TPC of the evolved populations were different from that of the control/collection population. (the two-curve model, fit1.5.P, had lower AICc than the one-curve model, fit0.5.P)


# 7/10/19 4:15 - updated results:

# [1] "5-1"
# dAICc df
# fit0.5.P  0.0  6 
# fit1.5.P  3.9  11

# [1] "5-2"
# dAICc df
# fit1.5.P  0.0  11
# fit0.5.P 16.9  6 

# [1] "5-3"
# dAICc df
# fit1.5.P  0.0  11
# fit0.5.P  0.8  6 

# [1] "5-4"
# dAICc df
# fit1.5.P  0.0  11
# fit0.5.P 19.7  6 


# remove legend
plot.5.list[[1]]<-plot.5.list[[1]]+theme(legend.position = 'none')
plot.5.list[[2]]<-plot.5.list[[2]]+theme(legend.position = 'none')
plot.5.list[[3]]<-plot.5.list[[3]]+theme(legend.position = 'none')
plot.5.list[[4]]<-plot.5.list[[4]]+theme(legend.position = 'none')

# export multi-panel figure
grid.arrange(plot.5.list[[1]],plot.5.list[[2]],plot.5.list[[3]],plot.5.list[[4]],nrow=2)
FS4<-arrangeGrob(plot.5.list[[1]],plot.5.list[[2]],plot.5.list[[3]],plot.5.list[[4]],nrow=2)

ggsave(paste("/Users/colin/Research/Active/Figure_S4.pdf",sep=''),FS4,width=8.25,height=7.75)



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Fig. S5: traits by evolution line at 100 Generations  ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# Use growthTools - package developed by CT Kremer
# for access to package, contact CT Kremer at kremerco@msu.edu
# package is to be publicly released on GitHub, ideally later in 2019
library(growthTools)

# Load growth rate data:
tmp100<-read.csv('growth_rates_at_100gen_042519.csv')

# split data sets for separate analyses
tmp100.L1<-tmp100 %>% filter(Evol.Strain %in% c('Collection','L1-1','L1-2','L1-3','L1-4'))
tmp100.5<-tmp100 %>% filter(Evol.Strain %in% c('Collection','5-1','5-2','5-3','5-4'))

# use growthTools to automate fitting of DE model TPC curves
fit.5<-tmp100.5 %>% group_by(Evol.Strain) %>% do(tpcs=get.decurve.tpc(.$Temperature,.$Growth.rate,plotQ = T,conf.bandQ = T,id=.$Evol.Strain))
fit.L1<-tmp100.L1 %>% group_by(Evol.Strain) %>% do(tpcs=get.decurve.tpc(.$Temperature,.$Growth.rate,plotQ = T,conf.bandQ = T,id=.$Evol.Strain))

# extract traits and confidence intervals from the DE curve fits:
fit.5.clean <- fit.5 %>% summarise(Evol.Strain,topt=tpcs$topt,b1=tpcs$b1,b2=tpcs$b2,d0=tpcs$d0,d2=tpcs$d2,s=tpcs$s,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,topt.lw=tpcs$ciF[1,1],topt.up=tpcs$ciF[1,2])
fit.L1.clean <- fit.L1 %>% summarise(Evol.Strain,topt=tpcs$topt,b1=tpcs$b1,b2=tpcs$b2,d0=tpcs$d0,d2=tpcs$d2,s=tpcs$s,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,topt.lw=tpcs$ciF[1,1],topt.up=tpcs$ciF[1,2])

# Calculate CI around umax for N-limited populations
umax.list<-rep(NA,nrow(fit.5.clean))
ci.list<-rep(NA,nrow(fit.5.clean))
for(i in 1:nrow(fit.5.clean)){
  # here's the umax:
  umax<-decurve(fit.5.clean$topt[i],fit.5.clean$topt[i],fit.5.clean$b1[i],fit.5.clean$b2[i],fit.5.clean$d0[i],fit.5.clean$d2[i])
  
  # figure out the confidence band around umax:
  xs<-fit.5.clean$topt[i]
  cfs<-c(fit.5.clean$topt[i],fit.5.clean$b1[i],fit.5.clean$b2[i],fit.5.clean$d0[i],fit.5.clean$d2[i],fit.5.clean$s[i])
  names(cfs)<-c('topt','b1','b2','d0','d2','s')
  dvs<-deltavar(fun=decurve(xs,topt,b1,b2,d0,d2),meanval=cfs,Sigma=fit.5[i,2][[1]][[1]]$vcov)
  ci.approx<-1.96*sqrt(dvs)
  ci.list[i]<-ci.approx
  umax.list[i]<-umax
}
fit.5.clean$umax<-umax.list
fit.5.clean$umax.lw<-umax.list-ci.list
fit.5.clean$umax.up<-umax.list+ci.list


# Calculate CI around umax for N-replete populations
umax.list<-rep(NA,nrow(fit.L1.clean))
ci.list<-rep(NA,nrow(fit.L1.clean))
for(i in 1:nrow(fit.L1.clean)){
  # here's the umax:
  umax<-decurve(fit.L1.clean$topt[i],fit.L1.clean$topt[i],fit.L1.clean$b1[i],fit.L1.clean$b2[i],fit.L1.clean$d0[i],fit.L1.clean$d2[i])
  
  # figure out the confidence band around umax:
  xs<-fit.L1.clean$topt[i]
  cfs<-c(fit.L1.clean$topt[i],fit.L1.clean$b1[i],fit.L1.clean$b2[i],fit.L1.clean$d0[i],fit.L1.clean$d2[i],fit.L1.clean$s[i])
  names(cfs)<-c('topt','b1','b2','d0','d2','s')
  dvs<-deltavar(fun=decurve(xs,topt,b1,b2,d0,d2),meanval=cfs,Sigma=fit.L1[i,2][[1]][[1]]$vcov)
  ci.approx<-1.96*sqrt(dvs)
  ci.list[i]<-ci.approx
  umax.list[i]<-umax
}
fit.L1.clean$umax<-umax.list
fit.L1.clean$umax.lw<-umax.list-ci.list
fit.L1.clean$umax.up<-umax.list+ci.list


# Generate plot based on the results:

f5<-unique(rbind(fit.5.clean,fit.L1.clean))
f5$group<-ifelse(grepl(f5$Evol.Strain,pattern = '5'),'N-limited evolved',ifelse(f5$Evol.Strain=='Collection','control','N-replete evolved'))

# boxes for confidence interval on Control population
f5 %>% filter(Evol.Strain=='Collection')
umax.lw<-0.7809177
umax.up<-0.8670009
topt.lw<-29.13706
topt.up<-30.70002
br1<-data.frame(topt=c(20,33,33,20),umax=c(umax.lw,umax.lw,umax.up,umax.up))
br2<-data.frame(topt=c(topt.lw,topt.up,topt.up,topt.lw),umax=c(0.5,0.5,1.2,1.2))

FS5<-f5 %>% filter(Evol.Strain!='Collection') %>%
  ggplot(aes(x=topt,y=umax))+
  geom_polygon(data=br1,fill=gray(0.9))+
  geom_polygon(data=br2,fill=gray(0.9))+
  geom_hline(yintercept = 0.824,linetype=2)+
  geom_vline(xintercept = 29.9,linetype=2)+
  geom_point(aes(group=Evol.Strain,colour=factor(group),shape=factor(group)),size=2.5)+
  geom_errorbar(aes(ymin=umax.lw,ymax=umax.up,colour=group),width=0.12)+
  geom_errorbarh(aes(xmin=topt.lw,xmax=topt.up,colour=group),height=0.01)+
  scale_color_manual('Treatment',values=c("#E69F00","#009E73"))+
  scale_shape_discrete('Treatment')+
  scale_x_continuous('Topt (C)')+
  scale_y_continuous('Maximum growth rate (1/day)')+
  coord_cartesian(xlim=c(27,32),ylim=c(0.65,1.16))+
  theme_bw()+
  theme(panel.grid = element_blank())
FS5

ggsave(filename = '/Users/colin/Research/Active/Figure_S5.pdf',plot = FS5,width=6,height=4.5)


# Also report Tmax estimates (no easy way to obtain CI's for Tmax):
f5[,c('Evol.Strain','tmax','tmin')]

#   Evol.Strain  tmax     tmin
#   5-1          33.7    7.50 
#   5-2          33.6    9.35 
#   5-3          33.7   -170.   
#   5-4          33.5   -1.46 
#   Collection   34.0    0.420
#   L1-1         33.8    5.72 
#   L1-2         33.7    2.31 
#   L1-3         34.1   -8.98 
#   L1-4         33.6    5.10


# Minor sanity check:

# check Topt from different sections are actually the same:
get.topt.from.model<-function(i,lst){
  sum(coef(lst[[i]])[c('topt','topt.T')])
}

# for N-limited 
sapply(seq(1,4,1),get.topt.from.model,f1s)
# 28.76666 29.22213 29.23782 28.38696

# for N-replete
sapply(seq(1,4,1),get.topt.from.model,f1s.L1)
#30.59741 29.62507 30.92367 29.15110

f5[,c('Evol.Strain','topt')]
# Evol.Strain  topt
# <fct>       <dbl>
# 1 5-1          28.8
# 2 5-2          29.2
# 3 5-3          29.4
# 4 5-4          28.4
# 5 Collection   29.9
# 6 L1-1         30.6
# 7 L1-2         29.6
# 8 L1-3         30.9
# 9 L1-4         29.1

# looks quite close



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### 200 Generation analysis ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# load phytoplankton data
tmp200<-read.xlsx2("200Generations_Maria_Data_curves.xlsx",1,colClasses=c("character","numeric","numeric"))

# decompose ID column into multiple treatment columns:
tmp200$ID<-as.character(tmp200$ID)
parts<-strsplit(tmp200$ID," +") 
tmp200$line<-unlist(lapply(parts, '[[', 1))
tmp200$temperature<-as.numeric(unlist(lapply(parts, '[[', 2)))
tmp200$Replicate<-unlist(lapply(parts, '[[', 3))
head(tmp200)

ilist200<-unique(tmp200$line) #A list of the unique IDs that need fitting



#### (1) Fit single DE model to all strains together: ####

topt.guess<-mean(tmp200$temperature[tmp200$Growth.rate==max(tmp200$Growth.rate)])

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,b1=NA,b2=NA,d0=NA,d2=NA,s=log(2))

# run fit:
fit0<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(temperature,topt,exp(b1),exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tmp200,control=list(maxit=10000))

# extract parameters for polished fit
cfg<-as.list(coef(fit0$res.best))
guesses<-as.list(cfg)
guesses<-list(topt=cfg$topt,b1=exp(cfg$b1),b2=exp(cfg$b2),d0=exp(cfg$d0),d2=exp(cfg$d2),s=exp(cfg$s))

# polish best fit model using formula interface
fit<-mle2(Growth.rate~dnorm(mean=decurve(temperature,topt,b1,b2,d0,d2),sd=s),
          start=guesses,data=tmp200,control=list(maxit=0))
summary(fit)

# predict curve from model fit
pd.df<-expand.grid(line=unique(tmp200$line),temperature=seq(10,35,0.1))
pd.df<-pd.df[order(pd.df$line,pd.df$temperature),]
pd.df$Growth.rate<-predict(fit,newdata=pd.df)
head(pd.df)

# visualize fits:
g1<-ggplot(tmp200,aes(x=temperature,y=Growth.rate))+
  geom_point(aes(shape=line))+
  geom_line(data=pd.df)+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage',breaks=c('ANCL1-1','Collect','L1-1-old','L1-3'),
                       labels=c('Ancestral','Control','replete N 1','replete N 2'))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  theme_bw()
g1

#ggsave(filename = '/Users/colin/Dropbox/EvolutionPaper_SharedFolder/Ecology Letters Review/modified analyses/200gen_fits_1curve_042519.pdf',plot = g1,width=6,height=4.5)


#### (2) Fit DE model by high N vs. controls: ####

tmp200$highN.Q<-ifelse(grepl(pattern = "L1-1-old",x = tmp200$line) | grepl(pattern = "L1-3",x = tmp200$line),1,0)
head(tmp200)

topt.guess<-mean(tmp200$temperature[tmp200$Growth.rate==max(tmp200$Growth.rate[tmp200$highN.Q==0]) & tmp200$highN.Q==0])

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,topt.T=0,b1=NA,b1.T=0,b2=NA,b2.T=0,d0=NA,d0.T=0,d2=NA,d2.T=0,s=log(2))

# run fit:
fit1<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(temperature,topt+topt.T*highN.Q,exp(b1+b1.T*highN.Q),exp(b2+b2.T*highN.Q),exp(d0+d0.T*highN.Q),exp(d2+d2.T*highN.Q)),sd=exp(s)),grids=grids,start=start,data=tmp200,control=list(maxit=10000))

# extract parameters for polished fit
cfg1<-as.list(coef(fit1$res.best))
guesses1<-list(topt=cfg1$topt,topt.T=cfg1$topt.T,
               b1=exp(cfg1$b1),b1.T=exp(cfg1$b1+cfg1$b1.T)-exp(cfg1$b1),
               b2=exp(cfg1$b2),b2.T=exp(cfg1$b2+cfg1$b2.T)-exp(cfg1$b2),
               d0=exp(cfg1$d0),d0.T=exp(cfg1$d0+cfg1$d0.T)-exp(cfg1$d0),
               d2=exp(cfg1$d2),d2.T=exp(cfg1$d2+cfg1$d2.T)-exp(cfg1$d2),
               s=exp(cfg1$s))

# polish best fit model using formula interface
fit1P<-mle2(Growth.rate~dnorm(mean=decurve(temperature,topt+topt.T*highN.Q,b1+b1.T*highN.Q,b2+b2.T*highN.Q,d0+d0.T*highN.Q,d2+d2.T*highN.Q),sd=s),
            start=guesses1,data=tmp200,control=list(maxit=0))
summary(fit1P)

# predict curves from fit
pd.df1<-merge(expand.grid(line=unique(tmp200$line),temperature=seq(10,35,0.1)),unique(tmp200[,c("line","highN.Q")]))
pd.df1<-pd.df1[order(pd.df1$line,pd.df1$temperature),]
pd.df1$Growth.rate<-predict(fit1P,newdata=pd.df1)
head(pd.df1)

# visualize fits:
g2<-ggplot(tmp200,aes(x=temperature,y=Growth.rate))+
  geom_point(aes(colour=factor(highN.Q),shape=line))+
  geom_line(data=pd.df1,aes(colour=factor(highN.Q)))+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage',breaks=c('ANCL1-1','Collect','L1-1-old','L1-3'),
                       labels=c('Ancestral','Control','replete N 1','replete N 2'))+
  scale_colour_manual('Treatment',values=c("black","#009E73"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  theme_bw()
g2

#ggsave(filename = '/Users/colin/Dropbox/EvolutionPaper_SharedFolder/Ecology Letters Review/modified analyses/200gen_fits_2curve_042519.pdf',plot = g2,width=6,height=4.5)



#### (3) Fit DE model by individual lineage: ####

# construct dummy variable columns
unique(tmp200$line)
tmp200$l1.Q<-ifelse(grepl(pattern = "L1-1-old",x = tmp200$line),1,0)
tmp200$l2.Q<-ifelse(grepl(pattern = "L1-3",x = tmp200$line),1,0)
tmp200$l3.Q<-ifelse(grepl(pattern = "Collect",x = tmp200$line),1,0)
head(tmp200)

topt.guess<-mean(tmp200$temperature[tmp200$Growth.rate==max(tmp200$Growth.rate[tmp200$line=='ANCL1-1']) & tmp200$line=='ANCL1-1'])

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,topt.l1=0,topt.l2=0,topt.l3=0,
            b1=NA,b1.l1=0,b1.l2=0,b1.l3=0,
            b2=NA,b2.l1=0,b2.l2=0,b2.l3=0,
            d0=NA,d0.l1=0,d0.l2=0,d0.l3=0,
            d2=NA,d2.l1=0,d2.l2=0,d2.l3=0,
            s=log(2))

# run fit:
fit2<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(temperature,
                                                         topt+topt.l1*l1.Q+topt.l2*l2.Q+topt.l3*l3.Q,
                                                         exp(b1+b1.l1*l1.Q+b1.l2*l2.Q+b1.l3*l3.Q),
                                                         exp(b2+b2.l1*l1.Q+b2.l2*l2.Q+b2.l3*l3.Q),
                                                         exp(d0+d0.l1*l1.Q+d0.l2*l2.Q+d0.l3*l3.Q),
                                                         exp(d2+d2.l1*l1.Q+d2.l2*l2.Q+d2.l3*l3.Q)),
                                            sd=exp(s)),grids=grids,start=start,data=tmp200,
                control=list(maxit=10000))

# extract parameters for polished fit
cfg2<-as.list(coef(fit2$res.best))
guesses2<-list(topt=cfg2$topt,topt.l1=cfg2$topt.l1,topt.l2=cfg2$topt.l2,topt.l3=cfg2$topt.l3,
               b1=exp(cfg2$b1),b1.l1=exp(cfg2$b1+cfg2$b1.l1)-exp(cfg2$b1),b1.l2=exp(cfg2$b1+cfg2$b1.l2)-exp(cfg2$b1),b1.l3=exp(cfg2$b1+cfg2$b1.l3)-exp(cfg2$b1),
               b2=exp(cfg2$b2),b2.l1=exp(cfg2$b2+cfg2$b2.l1)-exp(cfg2$b2),b2.l2=exp(cfg2$b2+cfg2$b2.l2)-exp(cfg2$b2),b2.l3=exp(cfg2$b2+cfg2$b2.l3)-exp(cfg2$b2),
               d0=exp(cfg2$d0),d0.l1=exp(cfg2$d0+cfg2$d0.l1)-exp(cfg2$d0),d0.l2=exp(cfg2$d0+cfg2$d0.l2)-exp(cfg2$d0),d0.l3=exp(cfg2$d0+cfg2$d0.l3)-exp(cfg2$d0),
               d2=exp(cfg2$d2),d2.l1=exp(cfg2$d2+cfg2$d2.l1)-exp(cfg2$d2),d2.l2=exp(cfg2$d2+cfg2$d2.l2)-exp(cfg2$d2),d2.l3=exp(cfg2$d2+cfg2$d2.l3)-exp(cfg2$d2),
               s=exp(cfg2$s))

# polish best fit model using formula interface
fit2P<-mle2(Growth.rate~dnorm(mean=decurve(temperature,topt+topt.l1*l1.Q+topt.l2*l2.Q+topt.l3*l3.Q,
                                           b1+b1.l1*l1.Q+b1.l2*l2.Q+b1.l3*l3.Q,
                                           b2+b2.l1*l1.Q+b2.l2*l2.Q+b2.l3*l3.Q,
                                           d0+d0.l1*l1.Q+d0.l2*l2.Q+d0.l3*l3.Q,
                                           d2+d2.l1*l1.Q+d2.l2*l2.Q+d2.l3*l3.Q),
                              sd=s),
            start=guesses2,data=tmp200,control=list(maxit=0))
summary(fit2P)

# extract predicted curves for plotting:
pd.df2<-merge(expand.grid(line=unique(tmp200$line),temperature=seq(10,35,0.1)),unique(tmp200[,c("line","highN.Q","l1.Q","l2.Q","l3.Q")]))
pd.df2<-pd.df2[order(pd.df2$line,pd.df2$temperature),]
pd.df2$Growth.rate<-predict(fit2P,newdata=pd.df2)
head(pd.df2)

# visualize fits:
g3<-ggplot(tmp200,aes(x=temperature,y=Growth.rate))+
  geom_point(aes(colour=factor(highN.Q),shape=line))+
  geom_line(data=pd.df2,aes(colour=factor(highN.Q),linetype=line))+
  #geom_line(data=pd.df1,aes(colour=factor(highN.Q)))+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Populations',breaks=c('ANCL1-1','Collect','L1-1-old','L1-3'),
                       labels=c('Ancestral','Control','N-Replete 1','N-Replete 2'))+
  scale_linetype_discrete('Populations',breaks=c('ANCL1-1','Collect','L1-1-old','L1-3'),
                          labels=c('Ancestral','Control','N-Replete 1','N-Replete 2'))+
  scale_colour_manual('Treatments',values=c("black","#009E73"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  theme_bw()
g3

#ggsave(filename = '/Users/colin/Dropbox/EvolutionPaper_SharedFolder/Ecology Letters Review/modified analyses/200gen_fits_4curve_042519.pdf',plot = g3,width=6,height=4.5)


# aggregate observations to avoid noisy plot
head(tmp200)

tmp200$Treatment<-ifelse(tmp200$line=='ANCL1-1','Ancestral',ifelse(tmp200$line=='Collect','Control','N-Replete'))
pd.df2$Treatment<-ifelse(pd.df2$line=='ANCL1-1','Ancestral',ifelse(pd.df2$line=='Collect','Control','N-Replete'))

  
tmp200.avg<-tmp200 %>% group_by(Treatment,line,highN.Q,temperature) %>% summarise(mn.Growth.rate=mean(Growth.rate),se.Growth.rate=ci.se(Growth.rate))

pd.df2 %>% group_by(line) %>% summarise(mx=max(Growth.rate))
# # A tibble: 4 x 2
#   line        mx
# 1 L1-1-old  1.16
# 2 L1-3      1.24
# 3 Collect   1.18
# 4 ANCL1-1   1.12


#### .......Create Fig. 1, panel c ####
F1.c<-ggplot(tmp200,aes(y=Growth.rate,x=temperature))+
  geom_point(data=tmp200.avg,aes(x=temperature,y=mn.Growth.rate,
                                 colour=Treatment,shape=Treatment))+
  geom_errorbar(data=tmp200.avg,aes(x=temperature,y=mn.Growth.rate,
                                   ymin=mn.Growth.rate-se.Growth.rate,
                                   ymax=mn.Growth.rate+se.Growth.rate,
                                   colour=Treatment),width=0.3,show.legend = F)+
  geom_line(data=pd.df2,aes(colour=Treatment,linetype=Treatment,group=line))+
  geom_hline(yintercept = 0)+
  scale_linetype_manual('',values=c(2,2,1))+
  scale_colour_manual('',values=c("gray","black","#009E73"))+
  scale_shape_manual('',values=c(4,4,1))+
  scale_x_continuous('Temperature (C)',limits=c(10,35))+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle('c)      ')



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### 200 generation model comparison and details ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

AICctab(fit,fit1P,fit2P,nobs=nrow(tmp200))
# dAICc df
# fit1P  0.0  11
# fit2P 22.1  21
# fit   93.9  6 


# allowing curves to differ between controls and high N lines improves fit to data
# but fitting individual curves to all 4 lines offers no improvement in fit

# comparison between 2-curve and 1-curve models is highly significant
anova(fit1P,fit)
# Tot Df Deviance  Chisq Df Pr(>Chisq)    
# 1     11  -88.787                         
# 2      6   17.341 106.13  5  < 2.2e-16 ***

# comparison between 4-curve and 1-curve models is not at all significant
anova(fit2P,fit1P)
# Tot Df Deviance  Chisq Df Pr(>Chisq)
# 1     21  -96.024                     
# 2     11  -88.787 7.2374 10     0.7029


# details of the best model:
summary(fit1P)

# Maximum likelihood estimation
# 
# Call:
#   mle2(minuslogl = Growth.rate ~ dnorm(mean = decurve(temperature, 
#                                                       topt + topt.T * highN.Q, b1 + b1.T * highN.Q, b2 + b2.T * 
#                                                         highN.Q, d0 + d0.T * highN.Q, d2 + d2.T * highN.Q), sd = s), 
#        start = guesses1, data = tmp200, control = list(maxit = 0))
# 
# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
#   topt   26.7776614  0.1326044 201.936 < 2.2e-16 ***
#   topt.T  4.2291911  0.0077907 542.850 < 2.2e-16 ***
#   b1      1.4280955  0.0403233  35.416 < 2.2e-16 ***
#   b1.T   28.5931720  0.0384039 744.538 < 2.2e-16 ***
#   b2      0.0918517  0.0018132  50.656 < 2.2e-16 ***
#   b2.T   -0.0898712  0.0018669 -48.139 < 2.2e-16 ***
#   d0      1.3879662  0.0771161  17.998 < 2.2e-16 ***
#   d0.T   29.2334607  0.0383844 761.596 < 2.2e-16 ***
#   d2      0.1083005  0.0012414  87.239 < 2.2e-16 ***
#   d2.T    0.5144806  0.0429789  11.970 < 2.2e-16 ***
#   s       0.1523699  0.0109952  13.858 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: -88.787 



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Fig. S6: 200 Generation TPCs, multi panel ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(cowplot)

g1b<-g1+theme(legend.position = 'none')
g2b<-g2+theme(legend.position = 'none')
g3b<-g3+theme(legend.position = 'none')
legend<-cowplot::get_legend(g3)

grid.arrange(g1b,g2b,g3b,legend,nrow=2)

FS6<-arrangeGrob(g1b,g2b,g3b,legend,nrow=2)
ggsave(filename = '/Users/colin/Research/Active/Figure_S6.pdf',plot = FS6,width=8.25,height=7.75)



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Fig. 1: TPCs a 100 and 200 generations ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Combine panels a-c into single figure with a hole in the bottom left:

# create fake data table for generating unified legend
dummy<-data.frame(Treatment=c('N-Limited Evolved','N-Replete Evolved','Ancestral','Control'),Growth.rate=1,Temperature=25)
dummy$Treatment<-factor(dummy$Treatment,levels=c('N-Limited Evolved','N-Replete Evolved','Ancestral','Control'))

leg.plot<-ggplot(dummy,aes(x=Temperature,y=Growth.rate))+
  geom_point(aes(colour=Treatment,shape=Treatment))+
  geom_line(aes(colour=Treatment,linetype=Treatment))+
  scale_linetype_manual('',values=c(1,1,2,2))+
  scale_colour_manual('',values=c("#E69F00","#009E73","gray","black"))+
  scale_shape_manual('',values=c(1,1,4,4))
legend.F1<-cowplot::get_legend(leg.plot)

# strip legends from individual panels
F1.a2<-F1.a+theme(legend.position = 'none')
F1.b2<-F1.b+theme(legend.position = 'none')
F1.c2<-F1.c+theme(legend.position = 'none')

# set up layout
hlay <- rbind(c(1,1,1,2,2,2),
              c(NA,3,NA,4,4,4))
grid.arrange(F1.a2,F1.b2,legend.F1,F1.c2,nrow=2,layout_matrix=hlay)

F1<-arrangeGrob(F1.a2,F1.b2,legend.F1,F1.c2,nrow=2,layout_matrix=hlay)
ggsave(filename = 'Figure_1.pdf',plot = F1,width=8.25,height=7.75)
