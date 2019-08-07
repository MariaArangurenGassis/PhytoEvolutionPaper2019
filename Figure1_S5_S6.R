
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

##### Estimate and compare the TPCs of evolved and control populations at 100 and 200 generations 

# Code by CT Kremer, last updated 4/25/19

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# load packages
library(xlsx)
library(lubridate)
library(dplyr)
library(bbmle)
library(ggplot2)
library(gridExtra)
library(devtools)
library(mleTools)

# Double exponential model, re-parameterized to depend explicitly on Topt:
decurve<-function(temp,topt,b1,b2,d0,d2){
  res <- b1*exp(b2*temp) - (d0 + ((b1*b2)/d2)*exp((b2-d2)*topt)*exp(d2*temp))
  res
}


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Calculate growth rates - 100 Generations ####

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# MARIA - this section is pretty complicated/messy for public release. I suggest generating a single .csv file with 100 generation growth rates, ideally using a separate R script, and then providing that as a single, already formatted input file to the remainder of this script.

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

# define a function for automating growth rate calculation:
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


write.xlsx2(rates100, file="100g_TPCs_GrowthRates_34_firstDays.xlsx",sheetName = "Sheet1",
            col.names = TRUE, row.names = TRUE, append = FALSE) #Added by Maria to extract the growth rates at 34ºC

#### Load previously calculated growth rates:

# load phytoplankton data
tmp100<-read.xlsx2("Curve100G_Chsim_NoPre.xlsx",1,colClasses=c("character","numeric","numeric","numeric","numeric"))
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

# take values from 34 C

# thin columns of tmp100
tmp100<-tmp100[,names(rates100)]

# take 34 C values from rates100 and combine
tmp100<-rbind(tmp100,rates100)
head(tmp100)

# save output:
write.csv(tmp100,'growth_rates_at_100gen_042519.csv',row.names=F)



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### 100 Generation analysis ####

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

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



#### ~~~~~~~~~~~~~ A. Low N results ~~~~~~~~~~~~~ ####

#### (1) Fit single DE model to all strains together: ####

topt.guess<-28.3

# set up a search grid of parameter guesses
grids<-list(b1=log(seq(0.01,5,0.5)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,5,0.5)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,b1=NA,b2=NA,d0=NA,d2=NA,s=log(5))

# run fit:
fit0<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt,exp(b1),exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tmp100.5,control=list(maxit=10000))

# extract parameters for polished fit
cfg<-as.list(coef(fit0$res.best))
guesses<-list(topt=cfg$topt,b1=exp(cfg$b1),b2=exp(cfg$b2),d0=exp(cfg$d0),d2=exp(cfg$d2),s=exp(cfg$s))

# polish best fit model using formula interface
fit.5<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
          start=guesses,data=tmp100.5,control=list(maxit=0))
summary(fit.5)

# predict curve shape, given model fit
pd.df.5<-expand.grid(Evol.Strain=unique(tmp100.5$Evol.Strain),Temperature=seq(10,34,0.1))
pd.df.5<-pd.df.5[order(pd.df.5$Evol.Strain,pd.df.5$Temperature),]
pd.df.5$Growth.rate<-predict(fit.5,newdata=pd.df.5)


# visualize fit:
p1.5<-ggplot(tmp100.5,aes(x=Temperature,y=Growth.rate))+
  geom_point(aes(shape=Evol.Strain))+
  geom_line(data=pd.df.5)+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage',breaks=c('Collection','5-1','5-2','5-3','5-4'),
                       labels=c('Control','5-1','5-2','5-3','5-4'))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1))+
  theme_bw()
p1.5

ggsave(filename = '100gen_lowN_fits_1curve_042519.pdf',plot = p1.5,width=6,height=4.5)




#### (2) Fit DE model by evolved vs. control: ####

topt.guess<-28

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,topt.T=0,b1=NA,b1.T=0,b2=NA,b2.T=0,d0=NA,d0.T=0,d2=NA,d2.T=0,s=log(2))

# run fit:
fit1<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,exp(b1+b1.T*evolved.Q),exp(b2+b2.T*evolved.Q),exp(d0+d0.T*evolved.Q),exp(d2+d2.T*evolved.Q)),sd=exp(s)),grids=grids,start=start,data=tmp100.5,control=list(maxit=10000))

# extract parameters for polished fit
cfg1<-as.list(coef(fit1$res.best))
guesses1<-list(topt=cfg1$topt,topt.T=cfg1$topt.T,
               b1=exp(cfg1$b1),b1.T=exp(cfg1$b1+cfg1$b1.T)-exp(cfg1$b1),
               b2=exp(cfg1$b2),b2.T=exp(cfg1$b2+cfg1$b2.T)-exp(cfg1$b2),
               d0=exp(cfg1$d0),d0.T=exp(cfg1$d0+cfg1$d0.T)-exp(cfg1$d0),
               d2=exp(cfg1$d2),d2.T=exp(cfg1$d2+cfg1$d2.T)-exp(cfg1$d2),
               s=exp(cfg1$s))

# polish best fit model using formula interface
fit1P.5<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,b1+b1.T*evolved.Q,b2+b2.T*evolved.Q,d0+d0.T*evolved.Q,d2+d2.T*evolved.Q),sd=s),
            start=guesses1,data=tmp100.5,control=list(maxit=0))
summary(fit1P.5)

# predict curves given model fit
pd.df1.5<-merge(expand.grid(Evol.Strain=unique(tmp100.5$Evol.Strain),Temperature=seq(10,34,0.1)),unique(tmp100.5[,c("Evol.Strain","evolved.Q")]))
pd.df1.5<-pd.df1.5[order(pd.df1.5$Evol.Strain,pd.df1.5$Temperature),]
pd.df1.5$Growth.rate<-predict(fit1P.5,newdata=pd.df1.5)
head(pd.df1.5)

# visualize fits:
p2.5<-ggplot(tmp100.5,aes(x=Temperature,y=Growth.rate))+
  geom_point(aes(colour=factor(evolved.Q),shape=Evol.Strain))+
  geom_line(data=pd.df1.5,aes(colour=factor(evolved.Q)))+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage')+
  scale_linetype_discrete('Lineage')+
  scale_colour_manual('Treatment',values=c("black","#E69F00"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1))+
  theme_bw()
p2.5

ggsave(filename = '100gen_lowN_fits_2curve_042519.pdf',plot = p2.5,width=6,height=4.5)


#### (3) Fit DE model by individual lineage: ####

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

# visualize fits:
p3.5<-ggplot(tmp100.5,aes(x=Temperature,y=Growth.rate))+
  geom_point(aes(colour=factor(evolved.Q),shape=Evol.Strain))+
  geom_line(data=pd.df2.5,aes(colour=factor(evolved.Q),linetype=Evol.Strain))+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage')+
  scale_linetype_discrete('Lineage')+
  scale_colour_manual('Treatment',values=c("black","#E69F00"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1))+
  theme_bw()
p3.5

ggsave(filename = '100gen_lowN_fits_5curve_042519.pdf',plot = p3.5,width=6,height=4.5)


#### Compare 100 generation low N models ####

# New interpretation: the 2 curve model is best, driven primarily by the behavior of the control at 12 C... sigh.
AICctab(fit.5,fit1P.5,fit2P.5,nobs=nrow(tmp100.5),weights=T)

# dAICc df weight
# fit1P.5  0.0  11 0.916 
# fit2P.5  5.2  26 0.069 
# fit.5    8.2  6  0.015 

# comparison between 2-curve and 1-curve models is significant
anova(fit1P.5,fit.5)
# Tot Df Deviance  Chisq Df Pr(>Chisq)   
# 1     11  -181.01                        
# 2      6  -161.36 19.646  5   0.001456 **

# comparison between 5-curve and 2-curve models is significant
anova(fit2P.5,fit1P.5)
# Tot Df Deviance  Chisq Df Pr(>Chisq)   
# 1     26  -216.20                        
# 2     11  -181.01 35.192 15   0.002309 **

# comparison between 5-curve and 1-curve models is highly significant
anova(fit2P.5,fit.5)
# Tot Df Deviance  Chisq Df Pr(>Chisq)    
# 1     26  -216.20                         
# 2      6  -161.36 54.838 20  4.341e-05 ***

# details of the best model:
summary(fit1P.5)

# Maximum likelihood estimation
# 
# Call:
#   mle2(minuslogl = Growth.rate ~ dnorm(mean = decurve(Temperature, 
#                                                       topt + topt.T * evolved.Q, b1 + b1.T * evolved.Q, b2 + b2.T * 
#                                                         evolved.Q, d0 + d0.T * evolved.Q, d2 + d2.T * evolved.Q), 
#                                        sd = s), start = guesses1, data = tmp100.5, control = list(maxit = 0))
# 
# Coefficients:
#   Estimate  Std. Error  z value     Pr(z)    
# topt   29.91543752  0.14206316 210.5784 < 2.2e-16 ***
#   topt.T -0.98041997  0.13748985  -7.1309 9.975e-13 ***
#   b1      3.78806555  0.04900981  77.2920 < 2.2e-16 ***
#   b1.T    6.04361457  0.05695678 106.1088 < 2.2e-16 ***
#   b2      0.00694467  0.00098321   7.0633 1.626e-12 ***
#   b2.T   -0.00320359  0.00100189  -3.1975 0.0013860 ** 
#   d0      3.79678027  0.04863515  78.0666 < 2.2e-16 ***
#   d0.T    6.19720139  0.05682720 109.0534 < 2.2e-16 ***
#   d2      0.77137238  0.05878886  13.1211 < 2.2e-16 ***
#   d2.T   -0.16733607  0.04537985  -3.6875 0.0002265 ***
#   s       0.12676622  0.00757565  16.7334 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: -181.0095





#### ~~~~~~~~~~~~~ B. High N results ~~~~~~~~~~~~~ ####

#### (1) Fit single DE model to all strains together: ####

topt.guess<-28.3

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,5,0.5)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,5,0.5)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,b1=NA,b2=NA,d0=NA,d2=NA,s=log(5))

# run fit:
fit0<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt,exp(b1),exp(b2),exp(d0),exp(d2)),sd=exp(s)),grids=grids,start=start,data=tmp100.L1,control=list(maxit=10000))

# extract parameters for polished fit
cfg<-as.list(coef(fit0$res.best))
guesses<-list(topt=cfg$topt,b1=exp(cfg$b1),b2=exp(cfg$b2),d0=exp(cfg$d0),d2=exp(cfg$d2),s=exp(cfg$s))

# polish best fit model using formula interface
fit.L1<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt,b1,b2,d0,d2),sd=s),
          start=guesses,data=tmp100.L1,control=list(maxit=0))
summary(fit.L1)

# predict curve from fit:
pd.df.L1<-expand.grid(Evol.Strain=unique(tmp100.L1$Evol.Strain),Temperature=seq(10,35,0.1))
pd.df.L1<-pd.df.L1[order(pd.df.L1$Evol.Strain,pd.df.L1$Temperature),]
pd.df.L1$Growth.rate<-predict(fit.L1,newdata=pd.df.L1)

# visualize fits:
p1.L1<-ggplot(tmp100.L1,aes(x=Temperature,y=Growth.rate))+
  geom_point(aes(shape=Evol.Strain))+
  geom_line(data=pd.df.L1)+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage',breaks=c('Collection','L1-1','L1-2','L1-3','L1-4'),
                       labels=c('Control','L1-1','L1-2','L1-3','L1-4'))+
  #scale_colour_manual('Treatment',values=c("black","#009E73"),
  #                    breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()
p1.L1

ggsave(filename = '100gen_highN_fits_1curve_042519.pdf',plot = p1.L1,width=6,height=4.5)


#### (2) Fit DE model by evolved vs. control: ####

topt.guess<-28

# set up search of a grid of parameter guesses
grids<-list(b1=log(seq(0.01,0.41,0.2)),b2=log(seq(0.1,0.5,0.2)),
            d0=log(seq(0.01,0.11,0.05)),d2=log(seq(0.1,0.7,0.2)))
start<-list(topt=topt.guess,topt.T=0,b1=NA,b1.T=0,b2=NA,b2.T=0,d0=NA,d0.T=0,d2=NA,d2.T=0,s=log(2))

# run fit:
fit1<-grid.mle2(minuslogl=Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,exp(b1+b1.T*evolved.Q),exp(b2+b2.T*evolved.Q),exp(d0+d0.T*evolved.Q),exp(d2+d2.T*evolved.Q)),sd=exp(s)),grids=grids,start=start,data=tmp100.L1,control=list(maxit=10000))

# extract parameters for polished fit
cfg1<-as.list(coef(fit1$res.best))
guesses1<-list(topt=cfg1$topt,topt.T=cfg1$topt.T,
               b1=exp(cfg1$b1),b1.T=exp(cfg1$b1+cfg1$b1.T)-exp(cfg1$b1),
               b2=exp(cfg1$b2),b2.T=exp(cfg1$b2+cfg1$b2.T)-exp(cfg1$b2),
               d0=exp(cfg1$d0),d0.T=exp(cfg1$d0+cfg1$d0.T)-exp(cfg1$d0),
               d2=exp(cfg1$d2),d2.T=exp(cfg1$d2+cfg1$d2.T)-exp(cfg1$d2),
               s=exp(cfg1$s))

# polish best fit model using formula interface
fit1P.L1<-mle2(Growth.rate~dnorm(mean=decurve(Temperature,topt+topt.T*evolved.Q,b1+b1.T*evolved.Q,b2+b2.T*evolved.Q,d0+d0.T*evolved.Q,d2+d2.T*evolved.Q),sd=s),
            start=guesses1,data=tmp100.L1,control=list(maxit=0))
summary(fit1P.L1)

# predict curves from model fit:
pd.df1.L1<-merge(expand.grid(Evol.Strain=unique(tmp100.L1$Evol.Strain),Temperature=seq(10,35,0.1)),unique(tmp100.L1[,c("Evol.Strain","evolved.Q")]))
pd.df1.L1<-pd.df1.L1[order(pd.df1.L1$Evol.Strain,pd.df1.L1$Temperature),]
pd.df1.L1$Growth.rate<-predict(fit1P.L1,newdata=pd.df1.L1)
head(pd.df1.L1)

# visualize fits:
p2.L1<-ggplot(tmp100.L1,aes(x=Temperature,y=Growth.rate))+
  geom_point(aes(colour=factor(evolved.Q),shape=Evol.Strain))+
  geom_line(data=pd.df1.L1,aes(colour=factor(evolved.Q)))+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage')+
  scale_colour_manual('Treatment',values=c("black","#009E73"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()
p2.L1

ggsave(filename = '100gen_highN_fits_2curve_042519.pdf',plot = p2.L1,width=6,height=4.5)



#### (3) Fit DE model by individual lineage: ####

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


# visualize fits:
p3.L1<-ggplot(tmp100.L1,aes(y=Growth.rate))+
  geom_point(aes(x=Temperature,colour=factor(evolved.Q),shape=Evol.Strain))+
  geom_line(data=pd.df2.L1,aes(x=Temperature,colour=factor(evolved.Q),linetype=Evol.Strain))+
  geom_hline(yintercept = 0)+
  scale_shape_discrete('Lineage')+  
  scale_linetype_discrete('Lineage')+
  scale_colour_manual('Treatment',values=c("black","#009E73"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()
p3.L1

ggsave(filename = '100gen_highN_fits_5curve_042519.pdf',plot = p3.L1,width=6,height=4.5)


#### Compare 100 generation high N models ####

# Now suggests that the 5-curve model is the best...
AICctab(fit.L1,fit1P.L1,fit2P.L1,nobs=nrow(tmp100.L1))
# dAICc df
# fit2P.L1  0.0  26
# fit.L1   51.0  6 
# fit1P.L1 54.6  11


# comparison between 2-curve and 1-curve models is not significant
anova(fit1P.L1,fit.L1)
# Tot Df Deviance  Chisq Df Pr(>Chisq)
# 1     11  -138.53                     
# 2      6  -130.77 7.7588  5       0.17

# comparison between 5-curve and 2-curve models is significant
anova(fit2P.L1,fit1P.L1)
#   Tot Df Deviance  Chisq Df Pr(>Chisq)    
# 1     26  -233.51                         
# 2     11  -138.53 94.985 15  1.155e-13 ***
  
# comparison between 5-curve and 1-curve models is significant
anova(fit2P.L1,fit.L1)
#   Tot Df Deviance  Chisq Df Pr(>Chisq)   
# 1     26  -233.51                         
# 2      6  -130.77 102.74 20  4.054e-13 ***

# details of the best model:
summary(fit2P.L1)

# Maximum likelihood estimation
# 
# Call:
#   mle2(minuslogl = Growth.rate ~ dnorm(mean = decurve(Temperature, 
#                                                       topt + topt.l1 * l1.Q + topt.l2 * l2.Q + topt.l3 * l3.Q + 
#                                                         topt.l4 * l4.Q, b1 + b1.l1 * l1.Q + b1.l2 * l2.Q + b1.l3 * 
#                                                         l3.Q + b1.l4 * l4.Q, b2 + b2.l1 * l1.Q + b2.l2 * l2.Q + 
#                                                         b2.l3 * l3.Q + b2.l4 * l4.Q, d0 + d0.l1 * l1.Q + d0.l2 * 
#                                                         l2.Q + d0.l3 * l3.Q + d0.l4 * l4.Q, d2 + d2.l1 * l1.Q + 
#                                                         d2.l2 * l2.Q + d2.l3 * l3.Q + d2.l4 * l4.Q), sd = s), 
#        start = guesses2, data = tmp100.L1, control = list(maxit = 0))
# 
# Coefficients:
#   Estimate  Std. Error  z value     Pr(z)    
# topt    29.91130370  0.13117256 228.0302 < 2.2e-16 ***
#   topt.l1  0.69553461  0.02930689  23.7328 < 2.2e-16 ***
#   topt.l2 -0.29666222  0.02016620 -14.7109 < 2.2e-16 ***
#   topt.l3  1.01146928  0.03232183  31.2937 < 2.2e-16 ***
#   topt.l4 -0.76089173  0.03147171 -24.1770 < 2.2e-16 ***
#   b1       4.71960929  0.04160863 113.4286 < 2.2e-16 ***
#   b1.l1    6.73353408  0.05702222 118.0861 < 2.2e-16 ***
#   b1.l2    2.23844271  0.05932336  37.7329 < 2.2e-16 ***
#   b1.l3    0.46704259  0.05674955   8.2299 < 2.2e-16 ***
#   b1.l4    2.98341608  0.06001698  49.7095 < 2.2e-16 ***
#   b2       0.00571558  0.00068300   8.3684 < 2.2e-16 ***
#   b2.l1   -0.00210639  0.00072125  -2.9205 0.0034950 ** 
#   b2.l2   -0.00165075  0.00083024  -1.9883 0.0467794 *  
#   b2.l3   -0.00202655  0.00088171  -2.2984 0.0215366 *  
#   b2.l4   -0.00028715  0.00080821  -0.3553 0.7223742    
#   d0       4.73463839  0.04116980 115.0027 < 2.2e-16 ***
#   d0.l1    6.95783992  0.05689145 122.3003 < 2.2e-16 ***
#   d0.l2    2.28911040  0.05918294  38.6786 < 2.2e-16 ***
#   d0.l3    0.28034677  0.05667308   4.9467 7.547e-07 ***
#   d0.l4    3.18400997  0.05976249  53.2777 < 2.2e-16 ***
#   d2       0.77336867  0.05354000  14.4447 < 2.2e-16 ***
#   d2.l1    0.26548525  0.04566577   5.8137 6.112e-09 ***
#   d2.l2   -0.01683338  0.04537574  -0.3710 0.7106543    
#   d2.l3    0.43388621  0.06277719   6.9115 4.795e-12 ***
#   d2.l4   -0.13678459  0.03922175  -3.4875 0.0004876 ***
#   s        0.10509317  0.00630476  16.6689 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: -233.5127




#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Fig. 1 a & b: 100 Generation TPC plots  ####

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# MARIA - this section of code uses the results of the 5-curve fits (above) to produce versions of Fig. 1a and 1b. I imagine you'll want to adjust the plot to further match the layout of your original Fig. 1, but hopefully this is a good start?

# calculate CI's based on estimated standard errors
ci.se<-function(x,na.rm=F){
  if(na.rm){
    x<-na.omit(x)
  }
  1.96*(sd(x)/length(x))
}

# To avoidmaking a cluttered plot, average across technical replicates within a population/temperature:
pnts.L1<-tmp100.L1 %>% group_by(Evol.Strain,evolved.Q,Temperature) %>% summarise(mn.gr=mean(Growth.rate),ci.gr=ci.se(Growth.rate))
str(pnts.L1)

fig1b<-ggplot(pnts.L1,aes(x=Temperature))+
  geom_line(data=pd.df2.L1,aes(y=Growth.rate,x=Temperature,colour=factor(evolved.Q),linetype=factor(evolved.Q),group=Evol.Strain),size=0.8)+
  geom_point(aes(y=mn.gr,colour=factor(evolved.Q),shape=factor(evolved.Q)),size=2.3)+
  geom_errorbar(aes(ymin=mn.gr-ci.gr,ymax=mn.gr+ci.gr,colour=factor(evolved.Q)),width=0.3)+
  geom_hline(yintercept = 0)+
  scale_shape_manual('Treatment',values=c(4,1),labels = c("Control \n(25 C)\n","N-replete evolved \n(31 C)"))+  
  scale_linetype_manual('Treatment',values=c(2,1),labels = c("Control \n(25 C)\n","N-replete evolved \n(31 C)"))+
  scale_colour_manual('Treatment',values=c("black","#009E73"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","N-replete evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()+
  theme(panel.grid = element_blank())
fig1b


# average across technical replicates within a population/temperature:
pnts.5<-tmp100.5 %>% group_by(Evol.Strain,evolved.Q,Temperature) %>% summarise(mn.gr=mean(Growth.rate),ci.gr=ci.se(Growth.rate))

fig1a<-ggplot(pnts.5,aes(x=Temperature))+
  geom_line(data=pd.df2.5,aes(y=Growth.rate,x=Temperature,colour=factor(evolved.Q),linetype=factor(evolved.Q),group=Evol.Strain),size=0.8)+
  geom_point(aes(y=mn.gr,colour=factor(evolved.Q),shape=factor(evolved.Q)),size=2.3)+
  geom_errorbar(aes(ymin=mn.gr-ci.gr,ymax=mn.gr+ci.gr,colour=factor(evolved.Q)),width=0.3)+
  geom_hline(yintercept = 0)+
  scale_shape_manual('Treatment',values=c(4,1),labels = c("Control \n(25 C)\n","N-limited evolved \n(31 C)"))+  
  scale_linetype_manual('Treatment',values=c(2,1),labels = c("Control \n(25 C)\n","N-limited evolved \n(31 C)"))+
  scale_colour_manual('Treatment',values=c("black","#E69F00"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","N-limited evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.5,1.3))+
  theme_bw()+
  theme(panel.grid = element_blank())
fig1a

# combined plot?
grid.arrange(fig1a,fig1b,nrow=1)

ggsave(filename = 'mock_Fig1a_042519.pdf',plot = fig1a,width=6,height=4.5)
ggsave(filename = 'mock_Fig1b_042519.pdf',plot = fig1b,width=6,height=4.5)



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Fig. 1 a & b: 100 Generation trait plots  ####

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

library(growthTools)
library(emdbook)

# use growthTools - package developed by CT Kremer
# for access to package, contact CT Kremer at kremerco@msu.edu
# package is to be released on GitHub, ideally later in 2019
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

suppFig<-f5 %>% filter(Evol.Strain!='Collection') %>%
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
suppFig

ggsave(filename = '100gen_umax_vs_topt_042519.pdf',plot = suppFig,width=6,height=4.5)


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



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### 200 Generation analysis ####

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

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

ggsave(filename = '200gen_fits_1curve_042519.pdf',plot = g1,width=6,height=4.5)


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

ggsave(filename = '200gen_fits_2curve_042519.pdf',plot = g2,width=6,height=4.5)



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
  scale_shape_discrete('Lineage',breaks=c('ANCL1-1','Collect','L1-1-old','L1-3'),
                       labels=c('Ancestral','Control','replete N 1','replete N 2'))+
  scale_linetype_discrete('Lineage',breaks=c('ANCL1-1','Collect','L1-1-old','L1-3'),
                       labels=c('Ancestral','Control','replete N 1','replete N 2'))+
  scale_colour_manual('Treatment',values=c("black","#009E73"),
                      breaks=c("0","1"),labels = c("Control \n(25 C)\n","Evolved \n(31 C)"))+
  scale_x_continuous('Temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  theme_bw()
g3

ggsave(filename = '200gen_fits_4curve_042519.pdf',plot = g3,width=6,height=4.5)




#### Compare 200 generation models ####

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


#######################################


