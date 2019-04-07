# replicate results from Jung et al (2014) and add some additional results
# this script deals mainly with endogeneity issues
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

# testing for endogeneity ####
# in their reply to Bakkensen&Larson Jung et al describe their test for endogeneity of NDAM as follows:
# "we built a simple original model in which fatality was regressed on normalized damage and gender index. 
# Next, we regressed normalized damage on gender index and minimum pressure and obtained residuals. 
# Finally, residuals were added as an additional regressor in the original model as well as in the count model. 
# The added residuals were not statistically different from zero."
model.endog.test.1=glm.nb(alldeaths~NDAM+MFI,data=data) 
summary(model.endog.test.1)
model.endog.test.2=lm(NDAM~MP+MFI,data=data) 
summary(model.endog.test.2)
data[,"res.endog.test.2"]=residuals(model.endog.test.2)
# now it's unclear what they mean by 'original model' and 'count model'. let's try two different possibilities
model.endog.test.3=glm.nb(alldeaths~NDAM+MFI+res.endog.test.2,data=data) 
summary(model.endog.test.3)
model.endog.test.4=glm.nb(alldeaths~MP+NDAM+MFI+MFI*MP+MFI*NDAM+res.endog.test.2,data=data) 
summary(model.endog.test.4)