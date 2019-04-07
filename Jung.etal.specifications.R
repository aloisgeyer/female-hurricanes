# replicate results from Jung et al (2014) and add some additional results
# this script deals with alternative specifications
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

# combine categories 4 and 5 (see below)
data[,"cat"]=data[,"Category"]
data[,"cat"]=4*(data[,"cat"]>3) + data[,"cat"]*(data[,"cat"]<4)

data[,"Category"]=factor(data[,"Category"])
data[,"cat"]=factor(data[,"cat"])

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

model.3=glm.nb(alldeaths~MP+NDAM+MFI+MFI*MP+MFI*NDAM,data=data)
summary(model.3)
lrtest(model.0,model.3)

model.4=glm.nb(alldeaths~ZMP+ZNDAM+ZMFI+ZMFI*ZMP+ZMFI*ZNDAM,data=data)
summary(model.4)
lrtest(model.0,model.4)
data[,"fit.4"]=fitted.values(model.4)

# alternative models ####
# LOGS of alldeaths ####
model.0.log=lm(log(1+alldeaths)~1,data=data)
model.1.log=lm(log(1+alldeaths)~MP,data=data)
summary(model.1.log)
lrtest(model.0.log,model.1.log)

model.2.log=lm(log(1+alldeaths)~MP+NDAM+MFI,data=data)
summary(model.2.log)
lrtest(model.0.log,model.2.log)

model.3.log=lm(log(1+alldeaths)~ZMP+ZNDAM+ZMFI+ZMFI*ZMP+ZMFI*ZNDAM,data=data)
summary(model.3.log)
lrtest(model.0.log,model.3.log)

# Gender_MF instead of MFI ####
model.3.dummy=glm.nb(alldeaths~MP+NDAM+Gender_MF+Gender_MF*MP+Gender_MF*NDAM,data=data)
summary(model.3.dummy)
lrtest(model.0,model.3.dummy)

model.4.dummy=glm.nb(alldeaths~ZMP+ZNDAM+Gender_MF+Gender_MF*ZMP+Gender_MF*ZNDAM,data=data)
summary(model.4.dummy)
lrtest(model.0,model.4.dummy)
data[,"fit.4.dummy"]=fitted.values(model.4.dummy)
summary(data$fit.4.dummy) # NOTE: two fitted values are STILL much larger than max in data! 
data %>% arrange(desc(fit.4.dummy)) %>% head() 

# EXCLUDE NDAM (endogenous?!) and replace it by 'Category' ####
# Category is not defined in the paper, but here: https://en.wikipedia.org/wiki/Saffir%E2%80%93Simpson_scale 

# descriptives ####
vars=c("alldeaths","MFI","MP","NDAM")
table(data[,"Category"])
data[,c(vars,"Category")] %>% group_by(Category) %>% summarize_at(.vars=vars,.funs=mean) %>% data.frame
data[,c(vars,"Category")] %>% group_by(Category) %>% summarize_at(.vars=vars,.funs=sd) %>% data.frame
data[,c(vars,"Category")] %>% group_by(Category) %>% summarize_at(.vars=vars,.funs=median) %>% data.frame

model.cat.1=glm.nb(alldeaths~MP+Category+MFI,data=data)
summary(model.cat.1)
lrtest(model.0,model.cat.1)
data[,"fit.cat.1"]=fitted.values(model.cat.1)

model.cat.2=glm.nb(alldeaths~MP+Category+MFI+MFI*MP+MFI*Category,data=data)
summary(model.cat.2)
lrtest(model.0,model.cat.2)
data[,"fit.cat.2"]=fitted.values(model.cat.2)

model.cat.3=glm.nb(alldeaths~MP+Category+MFI+MFI*MP,data=data)
summary(model.cat.3)
lrtest(model.0,model.cat.3)
data[,"fit.cat.3"]=fitted.values(model.cat.3)

# use reduced number of categories
model.cat.4=glm.nb(alldeaths~MP+cat+MFI+MFI*MP,data=data)
summary(model.cat.4)
lrtest(model.0,model.cat.4)
data[,"fit.cat.4"]=fitted.values(model.cat.4)

data.actual=cbind(data[,c("alldeaths","MFI")],"actual");names(data.actual)=c("alldeaths","MFI","type")
data.fitted=cbind(data[,c("fit.cat.3","MFI")],"fitted");names(data.fitted)=c("alldeaths","MFI","type")
data.4.plot=rbind(data.actual,data.fitted)
plt=ggplot(data.4.plot, aes(x=MFI,y=alldeaths,col=type) ) + geom_point()
plt = plt + scale_y_continuous(trans='log',limits=c(1,round(max(data$fit.4),-3))) 
plt

lo.qtl=0.2
hi.qtl=0.8

# reproduce Fig.1. using model.cat.3
data.fig1=NULL
data.fig1$MFI=seq(1,11,1)
data.fig1$MP=as.vector(quantile(data[,"MP"],lo.qtl))
data.fig1$Category=as.vector(1)
data.fig1=data.frame(data.fig1)
data.fig1=rbind(data.fig1,data.fig1,data.fig1,data.fig1,data.fig1)
data.fig1[12:22,"Category"]=as.vector(2)
data.fig1[23:33,"Category"]=as.vector(3)
data.fig1[34:44,"Category"]=as.vector(4)
data.fig1[45:55,"Category"]=as.vector(5)
data.fig1[,"Category"]=factor(data.fig1[,"Category"])
data.fig1[,"pred.fig1.cat.3.mp.lo"]=exp(predict(model.cat.3,newdata=data.fig1))

plt=ggplot(data.fig1, aes(x=MFI,y=pred.fig1.cat.3.mp.lo,col=Category) ) + geom_point()  + 
  labs(y="predicted fatality counts model.cat.3") + scale_y_continuous(trans='log',limits=c(1,150)) 
plt

# consider high level of MP
data.fig1[,"MP"]=as.vector(quantile(data[,"MP"],hi.qtl))
data.fig1[,"pred.fig1.cat.3.mp.hi"]=exp(predict(model.cat.3,newdata=data.fig1))

plt=ggplot(data.fig1, aes(x=MFI,y=pred.fig1.cat.3.mp.hi,col=Category) ) + geom_point()  + 
  labs(y="predicted fatality counts model.cat.3") + scale_y_continuous(trans='log',limits=c(1,150)) 
plt
