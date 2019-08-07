
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Code for determining model parameters given empirical growth rates ####

# By CT Kremer, 2019

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Load & format data ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(lubridate)
library(dplyr)

# load data from 200 generations
sat<-read.csv("/Users/colin/Dropbox/EvolutionPaper_SharedFolder/growth_rate_analyses/data/curves_200gen.csv")

# set up date column:
sat$date2<-mdy_hm(sat$date)

# naive calculation of delta time:
min.time<-min(sat$date2,na.rm = T)
sat$dtime<-difftime(sat$date2,min.time,units=c("days"))

# format
sat$evoln<-ifelse(grepl(pattern = "L1",sat$strain),"L1","Collection")
sat$temperature.f<-factor(sat$temperature)
sat$ids<-paste(sat$ids,sat$temperature)

# calculate 200 generation growth rates at each temperature for each population
uids<-unique(sat$ids)
rs<-rep(NA,length(uids))
for(i in 1:length(uids)){
  tmp<-sat[sat$ids==uids[i],]
  tmp<-tmp[order(tmp$dtime),]
  tmp<-tmp[1:5,]
  rs[i]<-coef(lm(log(fluor)~dtime,data=tmp))[[2]]
}
sat2<-merge(data.frame(ids=uids,r=rs),unique(sat[,c('ids','temperature','strain')]))
head(sat2)

# Visualize abundance time series supporting growth rate calculations
ggplot(sat,aes(x=dtime,y=log(fluor)))+
  geom_point(aes(colour=temperature.f,group=rep.flask))+
  geom_line(aes(colour=temperature.f,group=rep.flask))+
  facet_grid(temperature~strain)+
  ggtitle("at 200 G")

ggplot(sat[sat$temperature %in% c(34,35),],aes(x=dtime,y=log(fluor)))+
  geom_point(aes(colour=strain))+
  geom_line(aes(colour=strain,group=ids))+
  ggtitle("200G at 34 and 35 C")+
  facet_wrap(~temperature.f)

# visualize growth rates vs. temperature
sat2$category<-as.numeric(sat2$strain %in% c('ANCL1-1','Collect'))
ggplot(sat2,aes(x=temperature,y=r))+
  geom_point(aes(colour=factor(category)),alpha=0.3)+
  scale_colour_manual(values=c('black','red'))+
  theme_bw()

sat3<-sat2

# categorize strains
sat3$case<-ifelse(sat3$strain %in% c('ANCL1-1','Collect'),'lowP','highP')
sat3$case.code<-ifelse(sat3$case=='lowP',0,1)

# remove growth rates at 35 C for control populations - beyond the thermal niche, negative growth rate estimates are unreliable.
sat3<-sat3[!(sat3$temperature==35 & sat3$case=='lowP'),]

# save output for plotting in mathematica:
#write.csv(sat3,"/Users/colin/Dropbox/EvolutionPaper_SharedFolder/Data/data_for_tpc_fits_to_parameterize_model.csv",row.names=F)

# Exclude data at 10 C before model fitting (see appendix 1):
sat3<-sat3[sat3$temperature>15,]


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Estimate model parameters by fitting TPCs ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### Assume growth limited by r2 (see Appendix 1)

# growth equation:
intr.T.r2<-function(temp,b1,b2,d1,d2){
  cf<-((b1)*exp(b2*temp))-((d1)*exp(d2*temp))
  cf
} # Where d1 varies among strains

# Plot data and come up with initial guesses for parameters
plot(r~temperature,data=sat3,col=sat3$case.code+1,ylim=c(-0.6,2),xlim=c(10,36))
curve(intr.T.r2(x,0.04,0.21,0.022,0.23),-10,35,col='black',add=T)
curve(intr.T.r2(x,0.04,0.21,0.022-0.002,0.23),-10,35,col='red',lty=3,add=T)


# Define negative log likelihood calculator
nll.intr.T.r2<-function(b1,b2,d1,dB,d2,s){
  dat<-sat3
  #  w<-1/dat$sig
  w<-1
  mns<-intr.T.r2(dat$temperature,exp(b1),exp(b2),exp(d1+dB*dat$case.code),exp(d2))
  liks<-dnorm(dat$r,mean=mns,sd=exp(s),log=T)
  -sum(w*liks)
}

# Series of guesses tried while fussing with the optimization to obtain parameter estimates leading to global max log lik:
#guess<-list(b1=-1.355,b2=-2.823,d1=-12.95,dB=-0.701,d2=-0.901,s=-2.4)
#guess<-list(b1=-1.12,b2=-3.04,d1=-17.3,dB=-0.947,d2=-0.634,s=-2.4)
#guess<-list(b1=-1.189,b2=-3.0,d1=-17.3,dB=-0.935,d2=-0.634,s=-1.8)
#guess<-list(b1=-1.159,b2=-3.01,d1=-17.5,dB=-0.95,d2=-0.622,s=-1.8)
guess<-list(b1=-1.196,b2=-2.99,d1=-17.18,dB=-0.93,d2=-0.64,s=-2)

# Fit function using MLE and above guesses:
m2c<-mle2(nll.intr.T.r2,start=guess,control=list(maxit=5000,reltol=1e-12,ndeps=rep(1e-4,length(guess))),data=sat3)

# profile the fit, check surface shape, and extract parameter estimates
pf2c<-profile(m2c)
plot(pf2c)
cfs2<-c(exp(coef(m2c)[1:2]),coef(m2c)[3:4],exp(coef(m2c)[5]))
summary(m2c)

# visualize the fit
plot(r~temperature,data=sat3,col=sat3$case.code+1,ylim=c(-0.6,2),xlim=c(10,36),
     main='R2 limited exp. growth')
curve(intr.T.r2(x,cfs2[[1]],cfs2[[2]],exp(cfs2[[3]]),cfs2[[5]]),-10,35,col='black',add=T)
curve(intr.T.r2(x,cfs2[[1]],cfs2[[2]],exp(cfs2[[3]]+cfs2[[4]]),cfs2[[5]]),-10,35,col='red',add=T)

summary(m2c)
#       Estimate Std. Error z value     Pr(z)    
#  b1  -1.149704   0.096258 -11.944 < 2.2e-16 ***
#  b2  -3.025201   0.071610 -42.246 < 2.2e-16 ***
#  d1 -17.731242   0.297999 -59.501 < 2.2e-16 ***
#  dB  -0.959769   0.040257 -23.841 < 2.2e-16 ***
#  d2  -0.611235   0.015815 -38.649 < 2.2e-16 ***
#  s   -2.287018   0.080064 -28.565 < 2.2e-16 ***


# Note that:

# b1 in the summary = \alpha in the mathematica model
#exp(-1.891107)
exp(-1.149704) # 0.3167305

# d1 in the summary = \delta1 in the mathematica model
exp(-17.731242) # 1.992597e-08

# d1+dB in summary = \delta2 in mathematica model
exp(-17.731242+-0.959769) # 7.631274e-09

# b2 = b2 in model
exp(-3.025201) #0.04854806

# d2 = d2 in model
exp(-0.611235) #0.5426802

