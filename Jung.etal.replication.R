# replicate results from Jung et al (2014) and add some additional results
library(MASS)
library(lmtest)
library(dplyr)
library(ggplot2)

# NOTE: any suggestions which help to make this code more efficient are very welcome! 

# data ####
data=read.csv("Jung.etal.data.csv",sep=";")

# relabel
data[,"MP"]=data[,"MinPressure_before"]
data[,"MFI"]=data[,"MasFem"]

data[,"ZMFI"]=data[,"ZMasFem"]
data[,"ZMP"]=data[,"ZMinPressure_A"]
data[,"Gender_MF"]=factor(data[,"Gender_MF"],labels=c("male","female"))
data[,"Category"]=factor(data[,"Category"])

vars=c("alldeaths","MFI","MP","NDAM","ZMFI","ZMP","ZNDAM")

# descriptives ####
summary(data[,vars])
data %>% summarize_at(.vars=vars,.funs=mean)
data %>% summarize_at(.vars=vars,.funs=sd)
cor(data[,vars])
table(data[,"Gender_MF"])
data[,c(vars,"Gender_MF")] %>% group_by(Gender_MF) %>% summarize_at(.vars=vars,.funs=mean) %>% data.frame
data[,c(vars,"Gender_MF")] %>% group_by(Gender_MF) %>% summarize_at(.vars=vars,.funs=sd) %>% data.frame
data[,c(vars,"Gender_MF")] %>% group_by(Gender_MF) %>% summarize_at(.vars=vars,.funs=median) %>% data.frame

hist(data[,"alldeaths"])
hist(log(1+data[,"alldeaths"]))
ggplot(data, aes(alldeaths, col=Gender_MF)) + geom_density()
plot(data[,c("MFI","alldeaths")])
plot(data[,c("MFI","NDAM")])

data[,c("Name",vars)] %>% arrange(desc(alldeaths)) %>% head() # most severe in terms of alldeaths
data[,c("Name",vars)] %>% arrange(desc(NDAM)) %>% head() # most severe in terms of NDAM
data[,c("Name",vars)] %>% arrange(desc(Category)) %>% head() # most severe in terms of Category

# models in paper ####
model.0=glm.nb(alldeaths~1,data=data)
model.1=glm.nb(alldeaths~MP,data=data)
summary(model.1)
lrtest(model.0,model.1)

model.1.upd=glm.nb(alldeaths~Minpressure_Updated.2014,data=data) # check whether using 'Updated' makes a difference
summary(model.1.upd)
model.1.upd$deviance-model.1.upd$null.deviance
lrtest(model.0,model.1.upd)

model.2=glm.nb(alldeaths~MP+NDAM+MFI,data=data)
summary(model.2)
lrtest(model.0,model.2)
data[,"fit.2"]=fitted.values(model.2)
data %>% arrange(desc(fit.2)) %>% head() 

model.3=glm.nb(alldeaths~MP+NDAM+MFI+MFI*MP+MFI*NDAM,data=data)
summary(model.3)
lrtest(model.0,model.3)

model.4=glm.nb(alldeaths~ZMP+ZNDAM+ZMFI+ZMFI*ZMP+ZMFI*ZNDAM,data=data)
summary(model.4)
lrtest(model.0,model.4)
data[,"fit.4"]=fitted.values(model.4)
summary(data$fit.4) # NOTE: two fitted values are much larger than max in data! 
data %>% arrange(desc(fit.4)) %>% head() 

# standardize variables (double-check)
data.scaled=data %>% mutate_at(.vars=vars[!vars%in%"alldeaths"],.funs=funs(scale(.,center=T,scale=T)))
model.4.sc=glm.nb(alldeaths~MP+NDAM+MFI+MFI*MP+MFI*NDAM,data=data.scaled)
summary(model.4.sc) # should be the same as model.4 above
lrtest(model.0,model.4.sc)

# plot(log(1+data[,"alldeaths"]),log(1+data[,"alldeaths"]))
# points(log(1+data[,"alldeaths"]),log(1+data[,"fit"]),col="red")

# plot actual and fitted model 2
data.actual=cbind(data[,c("alldeaths","MFI")],"actual");names(data.actual)=c("alldeaths","MFI","type")
data.fitted=cbind(data[,c("fit.2","MFI")],"fitted");names(data.fitted)=c("alldeaths","MFI","type")
data.4.plot=rbind(data.actual,data.fitted)
plt=ggplot(data.4.plot, aes(x=MFI,y=alldeaths,col=type) ) + geom_point()
plt = plt + scale_y_continuous(trans='log',limits=c(1,round(max(data$fit.4),-3))) 
plt

# plot actual and fitted model 4
data.fitted=cbind(data[,c("fit.4","MFI")],"fitted");names(data.fitted)=c("alldeaths","MFI","type")
data.4.plot=rbind(data.actual,data.fitted)
plt=ggplot(data.4.plot, aes(x=MFI,y=alldeaths,col=type) ) + geom_point()
plt = plt + scale_y_continuous(trans='log',limits=c(1,round(max(data$fit.4),-3))) 
plt

# reproduce Fig.1 ####
# run two separate models ( see p.8786 top)
data[,"NDAM.cat"]=1*(data[,"NDAM"]>=quantile(data[,"NDAM"],0.5))
model.fig1=glm.nb(alldeaths~MP+NDAM.cat+MFI+MFI*MP+MFI*NDAM.cat,data=data, # does not converge without adding glm.control
                  control=glm.control(epsilon=1e-5,maxit=1000,trace=F))
summary(model.fig1)

data.fig1=NULL
data.fig1$MFI=seq(1,11,1)
data.fig1$ZMFI=(data.fig1$MFI-mean(data$MFI))/sd(data$MFI)
data.fig1$MP=mean(data$MP)
data.fig1$ZMP=0
data.fig1$NDAM.cat=0
data.fig1=data.frame(data.fig1)
data.fig1=rbind(data.fig1,data.fig1)
data.fig1[12:22,"NDAM.cat"]=1
data.fig1[,"NDAM"]=data.fig1[,"NDAM.cat"]
data.fig1[,"NDAM"]=factor(data.fig1[,"NDAM"],labels = c("below median","above median"))
data.fig1[,"pred.fig1"]=exp(predict(model.fig1,newdata=data.fig1))
max.pred.fig1=round(max(data.fig1$pred.fig1))

plt=ggplot(data.fig1, aes(x=MFI,y=pred.fig1,col=NDAM) ) + geom_point() + 
  labs(y="predicted fatality counts") + scale_y_continuous(limits=c(0,max.pred.fig1)) 
plt

# use original models for the same purpose; plug in lo.qtl%- and hi.qtl%-quantiles of NDAM 
# e.g. quantiles lo.qtl=0.25 and hi.qtl=0.75 can be viewed as the midpoints of the high-low category used in the paper
lo.qtl=0.2
hi.qtl=0.8
data.fig1=NULL
data.fig1$MFI=seq(1,11,1)
data.fig1$ZMFI=(data.fig1$MFI-mean(data$MFI))/sd(data$MFI)
data.fig1$MP=mean(data$MP)
data.fig1$ZMP=0
data.fig1$NDAM=as.vector(quantile(data[,"NDAM"],lo.qtl))
data.fig1=data.frame(data.fig1)
data.fig1=rbind(data.fig1,data.fig1)
data.fig1[12:22,"NDAM"]=as.vector(quantile(data[,"NDAM"],hi.qtl))
data.fig1$ZNDAM=(data.fig1$NDAM-mean(data$NDAM))/sd(data$NDAM)
data.fig1[,"pred.fig1.3"]=exp(predict(model.3,newdata=data.fig1)) # double check whether standardizing 
data.fig1[,"pred.fig1.4"]=exp(predict(model.4,newdata=data.fig1)) # makes a differnce  (it should not)
data.fig1[,"NDAM"]=factor(data.fig1[,"NDAM"],labels = c("NDAM low","NDAM high"))

plt=ggplot(data.fig1, aes(x=MFI,y=pred.fig1.3,col=NDAM) ) + geom_point() + 
  labs(y="predicted fatality counts from model 3") + scale_y_continuous(limits=c(0,max.pred.fig1)) 
plt

plt=ggplot(data.fig1, aes(x=MFI,y=pred.fig1.4,col=NDAM) ) + geom_point() + 
  labs(y="predicted fatality counts from model 4") + scale_y_continuous(limits=c(0,max.pred.fig1)) 
plt

# illustrate Christensen&Christensen's point regarding the need to account for both interaction terms
data.fig1=NULL
data.fig1$MFI=seq(1,11,1)
data.fig1$ZMFI=(data.fig1$MFI-mean(data$MFI))/sd(data$MFI)
data.fig1$MP=as.vector(quantile(data[,"MP"],lo.qtl))
data.fig1$NDAM=as.vector(quantile(data[,"NDAM"],lo.qtl))
data.fig1$type=1
data.fig1=data.frame(data.fig1)
data.fig1=rbind(data.fig1,data.fig1,data.fig1,data.fig1)
data.fig1[12:22,"MP"]=as.vector(quantile(data[,"MP"],hi.qtl))
data.fig1[12:22,"type"]=2
data.fig1[23:44,"NDAM"]=as.vector(quantile(data[,"NDAM"],hi.qtl))
data.fig1[23:44,"type"]=3
data.fig1[34:44,"MP"]=as.vector(quantile(data[,"MP"],hi.qtl))
data.fig1[34:44,"type"]=4
data.fig1$ZMP=(data.fig1$MP-mean(data$MP))/sd(data$MP)
data.fig1$ZNDAM=(data.fig1$NDAM-mean(data$NDAM))/sd(data$NDAM)
data.fig1[,"pred.fig1.3"]=exp(predict(model.3,newdata=data.fig1))
data.fig1[,"pred.fig1.4"]=exp(predict(model.4,newdata=data.fig1))
data.fig1[,"type"]=factor(data.fig1[,"type"],
                      labels = c("MP low; NDAM low","MP high; NDAM low","MP low; NDAM high","MP high; NDAM high"))

plt=ggplot(data.fig1, aes(x=MFI,y=pred.fig1.3,col=type) ) + geom_point()  + 
  labs(y="predicted fatality counts from model 3") + scale_y_continuous(limits=c(0,max.pred.fig1)) 
plt
plt=ggplot(data.fig1, aes(x=MFI,y=pred.fig1.4,col=type) ) + geom_point()  + 
  labs(y="predicted fatality counts from model 4") + scale_y_continuous(limits=c(0,max.pred.fig1)) 
plt

# additional models ####
# Malter claims that an additional interaction term should be included
malter.2=glm.nb(alldeaths~MP+NDAM+MP*NDAM+MFI+MFI*MP,data=data, # does not converge without adding glm.control
                control=glm.control(epsilon=1e-5,maxit=1000,trace=F)) # does not converge for epsilon>1.e-5 or maxit too low
summary(malter.2)
lrtest(model.0,malter.2)

malter.3=glm.nb(alldeaths~MP+NDAM+MP*NDAM+MFI+MFI*NDAM,data=data, # default setting does not converge
                control=glm.control(epsilon=1e-5,maxit=1000,trace=F))
summary(malter.3)
lrtest(model.0,malter.3)

malter.4=glm.nb(alldeaths~MP+NDAM+MP*NDAM+MFI+MFI*MP+MFI*NDAM,data=data, 
                control=glm.control(epsilon=1e-3,maxit=1000,trace=F)) # CANNOT FIND a SETTING that CONVERGES! 
# summary(malter.4)
# lrtest(model.0,malter.4)
