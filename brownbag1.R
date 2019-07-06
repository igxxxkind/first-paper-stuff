#install.packages("mFilter")
library(mFilter)

library(xtable)
library(fpp)
library(readxl)
library(ggplot2)
library(vars)
library(seasonal)
library(readxl)
library(graphics)
library(lmtest)
library(MARX)

setwd("D://Data//PhD//BrownBag")
source("functions.R")
unrate = data.frame(read_excel("UNRATE_M1953-2019.xls",sheet = "data"))[,2]

adf.test(unrate) #reject nonstationarity at 5%
pp.test(unrate)  #reject stationarity
pp.test(unrate,type = "Z(t_alpha)") #reject stationarity
####HP-filter####

hpfiltered = hpfilter(x = unrate,type = "lambda",freq = 14400)
#hpfiltered2 = hpfilter(x = unrate,type = "lambda",freq = 28800)


ts.plot(unrate)
lines(hpfiltered$trend, col = "red")
lines(hpfiltered$cycle+6, col = "blue")
abline(a=6,b = 0)
var(hpfiltered$cycle)

# ts.plot(unrate)
# lines(hpfiltered2$trend, col = "red")
# lines(hpfiltered2$cycle+6, col = "blue")
# abline(a=6,b = 0)
# var(hpfiltered2$cycle)
adf.test(hpfiltered$cycle)

###Baxter-King filter####

bkfiltered1 = bkfilter(unrate,pl = 2,pu = 96,nfix = 36) 
#filter out the irregular component and the business cycle component

bkfiltered2 = bkfilter(unrate,pl = 18,pu = 96,nfix = 36) 
#filter out the business cycle of 1.5-8 years

bkfiltered11 = bkfilter(unrate,pl = 2,pu = 96,nfix = 24) 
#filter out the business cycle of 1.5-8 years

bkfiltered12 = bkfilter(unrate,pl = 18,pu = 96,nfix = 24) 
#filter out the business cycle of 1.5-8 years


ts.plot(bkfiltered1$cycle, col = "red")
lines(bkfiltered11$cycle, col = "blue")
#lines(bkfiltered$cycle, col = "green4")
abline(0,0)

ts.plot(bkfiltered2$cycle, col = "red")
lines(bkfiltered12$cycle, col = "blue")
#lines(bkfiltered$cycle, col = "green4")
abline(0,0)

ts.plot(bkfiltered2$trend, col = "red")
lines(bkfiltered12$trend, col = "blue")
#lines(bkfiltered$cycle, col = "green4")
abline(0,0)

###Filtration is done###
####Step 2. Data splitting and lag selection####
unrate.ts = ts(data = unrate,start = c(1953,1),end = c(2019,5),frequency = 12)

cpi = data.frame(read_excel("CPIAUCSL_M1952-2019.xls",sheet = "data"))[,2]
cpi_log = ts(data = diff(log(cpi)),start = c(1952,1),end = c(2019,4),frequency = 12)
cpi.ts = window(x = cpi_log,start = c(1953,1),end = c(2019,4))
cpi.gr = ts(data = cpi[-1]/cpi[-length(cpi)],start = c(1952,1),end = c(2019,4),frequency = 12)
cpi.grate = window(x = cpi.gr,start = c(1953,1),end = c(2019,4))

unratecycle.ts = ts(data = bkfiltered12$cycle,start = c(1953,1),end = c(2019,5),frequency = 12)
unratecycle2.ts = ts(data = bkfiltered11$cycle,start = c(1953,1),end = c(2019,5),frequency = 12)

cpi1.ts = cpi.ts[!is.na(unratecycle.ts[-length(unratecycle.ts)])]
cpi1.gr = cpi.grate[!is.na(unratecycle.ts[-length(unratecycle.ts)])]

unratecycle.ts = unratecycle.ts[!is.na(unratecycle.ts)]
unratecycle2.ts = unratecycle2.ts[!is.na(unratecycle2.ts)]

ts.plot(cpi.ts*100)
lines(unratecycle.ts, col = "red")
lines(unratecycle2.ts, col = "blue")


unratecycle.ts = ts(data = unratecycle.ts,start = c(1955,1),end = c(2017,5),frequency = 12)
unratecycle2.ts = ts(data = unratecycle2.ts,start = c(1955,1),end = c(2017,5),frequency = 12)

cpi2.ts = ts(data = cpi1.ts,start = c(1955,1),end = c(2017,5),frequency = 12)
cpi2.gr = ts(data = cpi1.gr,start = c(1955,1),end = c(2017,5),frequency = 12)


###Divide the data into three major periods
##Great Inflation: 1960 - 1983
GGI = list()
GGI$cpi = window(x = cpi2.ts, start = c(1960,1), end = c(1982,12))
GGI$cpi2= window(x = cpi2.ts, start = c(1965,1), end = c(1982,12))
GGI$cpi = GGI$cpi2


GGI$unrate1 = window(x = unratecycle.ts, start = c(1960,1), end = c(1982,12))
GGI$unrate3 = window(x = unratecycle.ts, start = c(1965,1), end = c(1982,12))

GGI$unrate = GGI$unrate3
GGI$n = length(GGI$cpi)
GGI$stdcpi = (GGI$cpi - mean(GGI$cpi))/sqrt(var(GGI$cpi))


#check = modelfit2(GGI$cpi,l = 20)
GGI$lagselect = modelfit2(GGI$cpi,l = 20,
                         exo = GGI$unrate)

GGI$unitroot0 = summary(ur.df(GGI$cpi,type = "none"))
GGI$unitrootD = summary(ur.df(GGI$cpi,type = "drift"))
GGI$unitrootT = summary(ur.df(GGI$cpi,type = "trend"))

GGI$causalmodel = arima(x = GGI$cpi,order = c(2,0,0),
                       xreg = GGI$unrate,method = "CSS")


GGI$ACtest1 = diagnostics(cpi = GGI$cpi, unrate = GGI$unrate,range = 20,order = 30)
GGI$ACtest2 = diagnostics(cpi = GGI$cpi, unrate = GGI$unrate,range = 20,order = 10)

ts.plot(GGI$ACtest1$JBp.value)

##Great Moderation: 1985 - 2007
GGM = list()
GGM$cpi = window(x = cpi2.ts, start = c(1983,1), end = c(2006,12))
GGM$unrate1 = window(x = unratecycle.ts, start = c(1983,1), end = c(2006,12))
GGM$unrate2 = window(x = unratecycle2.ts, start = c(1983,1), end = c(2006,12))
GGM$n = length(GGM$cpi)
GGM$stdcpi = (GGM$cpi - mean(GGM$cpi))/sqrt(var(GGM$cpi))


GGM$lagselect = modelfit2(GGM$cpi,l = 20,exo = GGM$unrate1)
GGM$lagselect2 = modelfit2(GGM$cpi,l = 20,exo = GGM$unrate2)

GGM$ACtest1 = diagnostics(cpi = GGM$cpi, unrate = GGM$unrate1,range = 20,order = 30)
GGM$ACtest2 = diagnostics(cpi = GGM$cpi, unrate = GGM$unrate1,range = 20,order = 10)

GGM$unitroot0 = summary(ur.df(GGM$cpi,type = "none"))
GGM$unitrootD = summary(ur.df(GGM$cpi,type = "drift"))
GGM$unitrootT = summary(ur.df(GGM$cpi,type = "trend"))



##Great Recession: 2007 - 2017
GGR = list()
GGR$cpi = window(x = cpi2.ts, start = c(2007,1), end = c(2017,5))

GGR$unrate1 = window(x = unratecycle.ts, start = c(2007,1), end = c(2017,5))
GGR$unrate2 = window(x = unratecycle2.ts, start = c(2007,1), end = c(2017,5))

GGR$n = length(GGR$cpi)

GGR$stdcpi = (GGR$cpi - mean(GGR$cpi))/sqrt(var(GGR$cpi))

GGR$lagselect = modelfit2(GGR$cpi,l = 20,exo = GGR$unrate1)
GGR$lagselect2 = modelfit2(GGR$cpi,l = 20,exo = GGR$unrate2)

GGR$ACtest1 = diagnostics(cpi = GGR$cpi, unrate = GGR$unrate1,range = 20,order = 30)
GGR$ACtest2 = diagnostics(cpi = GGR$cpi, unrate = GGR$unrate1,range = 20,order = 10)

GGR$unitroot0 = summary(ur.df(GGR$cpi,type = "none"))
GGR$unitrootD = summary(ur.df(GGR$cpi,type = "drift"))
GGR$unitrootT = summary(ur.df(GGR$cpi,type = "trend"))

####Tables####
library(e1071)  
tables = list()

tables$sumstat = matrix(0,10,3)
tables$sumstat[,1] = c(mean(GGI$cpi),median(GGI$cpi),var(GGI$cpi),skewness(GGI$cpi),
                       kurtosis(GGI$cpi),min(GGI$cpi),max(GGI$cpi),
                       GGI$unitroot0@teststat,GGI$unitrootD@teststat[1],
                       GGI$unitrootT@teststat[1])
tables$sumstat[,2] = c(mean(GGM$cpi),median(GGM$cpi),var(GGM$cpi),skewness(GGM$cpi),
                       kurtosis(GGM$cpi),min(GGM$cpi),max(GGM$cpi),
                       GGM$unitroot0@teststat,GGM$unitrootD@teststat[1],
                       GGM$unitrootT@teststat[1])
tables$sumstat[,3] = c(mean(GGR$cpi),median(GGR$cpi),var(GGR$cpi),skewness(GGR$cpi),
                       kurtosis(GGR$cpi),min(GGR$cpi),max(GGR$cpi),
                       GGR$unitroot0@teststat,GGR$unitrootD@teststat[1],
                       GGR$unitrootT@teststat[1])
rownames(tables$sumstat) = c("Mean","Median","Variance","Skewness","Kurtosis",
                             "Min", "Max", "ADF test (origin)", 
                             "ADF test (const)","ADF test (const+trend)")

tables$tests = matrix(0,18,21)
tables$tests[1:6,] = rbind(GGI$ACtest1$LBtest,GGI$ACtest1$LBp.value,
                     GGI$ACtest1$BGtest,GGI$ACtest1$BGp.value,
                     GGI$ACtest1$JBstat,GGI$ACtest1$JBp.value)
tables$tests[7:12,] = rbind(GGM$ACtest1$LBtest,GGM$ACtest1$LBp.value,
                     GGM$ACtest1$BGtest,GGM$ACtest1$BGp.value,
                     GGM$ACtest1$JBstat,GGM$ACtest1$JBp.value)
tables$tests[13:18,] = rbind(GGR$ACtest1$LBtest,GGR$ACtest1$LBp.value,
                     GGR$ACtest1$BGtest,GGR$ACtest1$BGp.value,
                     GGR$ACtest1$JBstat,GGR$ACtest1$JBp.value)
tables$testrownames = c("Ljung-Box, ord=30", "L-B p.value",
                        "Breusch-Godfrey, ord=30", "B-G p.value",
                        "Jarque-Bera", "J-B p.value")
rownames(tables$tests) = rep(tables$testrownames,3)
colnames(tables$tests) = c(1:21)

tables$tests2 = matrix(0,18,21)
tables$tests2[1:6,] = rbind(GGI$ACtest2$LBtest,GGI$ACtest2$LBp.value,
                           GGI$ACtest2$BGtest,GGI$ACtest2$BGp.value,
                           GGI$ACtest2$JBstat,GGI$ACtest2$JBp.value)
tables$tests2[7:12,] = rbind(GGM$ACtest2$LBtest,GGM$ACtest2$LBp.value,
                            GGM$ACtest2$BGtest,GGM$ACtest2$BGp.value,
                            GGM$ACtest2$JBstat,GGM$ACtest2$JBp.value)
tables$tests2[13:18,] = rbind(GGR$ACtest2$LBtest,GGR$ACtest2$LBp.value,
                             GGR$ACtest2$BGtest,GGR$ACtest2$BGp.value,
                             GGR$ACtest2$JBstat,GGR$ACtest2$JBp.value)
tables$testrownames = c("Ljung-Box, ord=10", "L-B p.value",
                        "Breusch-Godfrey, ord=10", "B-G p.value",
                        "Jarque-Bera", "J-B p.value")
rownames(tables$tests2) = rep(tables$testrownames,3)
colnames(tables$tests2) = c(1:21)


tables$aicbic = matrix(0,6,20)
tables$aicbic[1:2,] = rbind(GGI$lagselect$aic,GGI$lagselect$bic)
tables$aicbic[3:4,] = rbind(GGM$lagselect$aic,GGM$lagselect$bic)
tables$aicbic[5:6,] = rbind(GGR$lagselect$aic,GGR$lagselect$bic)
rownames(tables$aicbic) = rep(c("AIC","BIC"),3)

####Step 3. Some tests####


##Breusch-Godfrey test
###works asymptotically well
y = GGI$cpi
n = length(y)

regs = matrix(GGI$unrate1,,1)
p=2
for (i in 1:p){
  #t = dim(regs)[1]
  regs = cbind(regs[-1,],y[-c((n-i+1):n)])
  }

ols = summary(lm(y[-c(1:p)]~1 + regs))



ord = 10
uhat = ols$residuals
m = length(uhat)
regs2 = regs

for (i in 1:ord){
  regs2 = cbind(regs2[-1,],uhat[-c((m-i+1):m)])
}

bgols = summary(lm(uhat[-c(1:ord)]~1 + regs2))
bgn = length(bgols$residuals)
#archols = 

R22 = 1 - sum(bgols$residuals^2)/sum((uhat[-c(1:ord)])^2)



LM = (bgn)*R22
FLM = R22/(1-R22) * (bgn-p-ord-1-1)/ord

bgtest(formula = y[-c(1:p)]~1 + regs,order = ord,type = "Chisq",fill = NA)
bgtest(formula = y[-c(1:p)]~1 + regs,order = ord,type = "F", fill = NA)

###ARCh-type test on heteroscedasticity
dim(regs)
uhatsq = (regs2[,-c(1:dim(regs)[2])])^2

archols = summary(lm(uhat[-c(1:ord)]~1 + uhatsq))

R2b=archols$r.squared
LMarch = R2b*length(archols$residuals)
pchisq(LMarch,df = ord)

###RESET
resettest(y[-c(1:p)]~1 + regs)

olsfv = (lm(y[-c(1:p)]~1 + regs))
yhat = fitted.values(olsfv)
yhat2 = cbind(yhat^2,yhat^3)
regs3 = cbind(regs,yhat2)

olsreset = lm(uhat~1+regs3) 
vhat = olsreset$residuals

resetstat = ((sum(uhat^2) - sum(vhat^2))/2)/(sum(vhat^2)/(length(uhat)-dim(regs)[2]+1-3-1))


#check = diagnostics(cpi = GGI$cpi, unrate = GGI$unrate1,range = 2,order = 30)

####Step 3. Analysis of the noncausal setings####
###Follow Lanne, Saikonen lag order selection

###Great Inflation Period###
n = length(GGI$cpi)
nlag = 12
GGI$lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                             GGI$cpi[1:(n - x)]))
#lags = lags[-c(1:nlag),]

GGI$AIC = c()
GGI$AIC2 = c()

for(i in 1:nlag){
exo = GGI$unrate1[-c(1:i)]
auxm = lm(GGI$cpi[-c(1:i)]~1+GGI$lags[-c(1:i),1:i]+exo)
GGI$AIC[i] = AIC(auxm)}

GGI$lags2 = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        GGI$stdcpi[1:(n - x)]))
for(i in 1:nlag){
  exo = GGI$unrate2[-c(1:i)]
  auxm = lm(GGI$stdcpi[-c(1:i)]~1+GGI$lags2[-c(1:i),1:i]+exo)
  GGI$AIC2[i] = AIC(auxm)}

#Lag selection process is a sequential procedure
#therefore simple commands from MARX package are not appropriate

#=>three possible candidates c(2,5,9)

GGI$lagorder = c(2,5,9)
exo2 = GGI$unrate1[-c(1:GGI$lagorder[1])]
GGI$cm2 = lm(GGI$cpi[-c(1:GGI$lagorder[1])]~1+
                        GGI$lags[-c(1:GGI$lagorder[1]),1:GGI$lagorder[1]]+exo2)


exo5 = GGI$unrate1[-c(1:GGI$lagorder[2])]
GGI$cm5 = lm(GGI$cpi[-c(1:GGI$lagorder[2])]~1+
                        GGI$lags[-c(1:GGI$lagorder[2]),1:GGI$lagorder[2]]+exo5)

exo9 = GGI$unrate1[-c(1:GGI$lagorder[3])]
GGI$cm9 = lm(GGI$cpi[-c(1:GGI$lagorder[3])]~1+
                        GGI$lags[-c(1:GGI$lagorder[3]),1:GGI$lagorder[3]]+exo9)


GGI$cmstat2 = summary(GGI$cm2)
GGI$cmstat5 = summary(GGI$cm5)
GGI$cmstat9 = summary(GGI$cm9)

GGI$MAR = list()
GGI$MAR$causal2 = selection.lag.lead(y = GGI$cpi, x = GGI$unrate,p_pseudo = 2)

GGI$MAR$causal5 = selection.lag.lead(y = GGI$cpi,x = GGI$unrate,p_pseudo = 5)

GGI$MAR$causal9 = selection.lag.lead(y = GGI$cpi,x = GGI$unrate,p_pseudo = 9)

GGI$MAR$model11 = mixed(y = GGI$cpi,x = GGI$unrate,p_C = 2,p_NC = 0)
GGI$MAR$model12 = mixed(y = GGI$cpi,x = GGI$unrate,p_C = 1,p_NC = 1)
GGI$MAR$model13 = mixed(y = GGI$cpi,x = GGI$unrate,p_C = 0,p_NC = 2)

GGI$MAR$table = matrix(0,6,8)
GGI$MAR$table[1,]= c(GGI$MAR$model11$coefficients[c(1,2,3,4,4,5,6,7)])
GGI$MAR$table[2,]= c(GGI$MAR$model11$se[c(1,2,3,4,4,5,6,7)])
GGI$MAR$table[3,]= c(GGI$MAR$model12$coefficients[c(1,2,2,3,3,4,5,6)])
GGI$MAR$table[4,]= c(GGI$MAR$model12$se[c(1,2,2,3,3,4,5,6)])
GGI$MAR$table[5,]= c(GGI$MAR$model13$coefficients[c(1,2,2,3,4,5,6,7)])
GGI$MAR$table[6,]= c(GGI$MAR$model13$se[c(1,2,2,3,4,5,6,7)])
rownames(GGI$MAR$table) = c("AR(2,0)","s.e.", "MAR(1,1)","s.e.", "AR(0,2)","s.e.")
colnames(GGI$MAR$table) = c("int", "lag1", "lag2", "lead1", "lead2", "exo", "df", "scale" )

GGI$MAR$model11$MSE = sum((GGI$cpi[-c(1,2)]-GGI$cm2$fitted.values)^2)
GGI$MAR$model12$MSE = sum((GGI$cpi[-c(1,n)]-GGI$MAR$model12$fitted.values)^2)
GGI$MAR$model13$MSE = sum((GGI$cpi[-c(n-1,n)]-GGI$MAR$model13$fitted.values)^2)

nlead = 2
GGI$leads = sapply(1:nlead, function(x) c(GGI$cpi[-c(1:x)],rep(NA, length.out = x)))


Box.test(GGI$MAR$model11$residuals,lag = 12,type = "Ljung-Box")
Box.test(GGI$MAR$model12$residuals,lag = 12,type = "Ljung-Box")
Box.test(GGI$MAR$model13$residuals,lag = 12,type = "Ljung-Box")

selection.lag.lead(GGI$cpi,x = GGI$unrate1,p_pseudo = 2)
GGI$MAR$select = selection.lag.lead(GGI$cpi2,x = GGI$unrate,p_pseudo = 2)

selection.lag.lead(GGI$stdcpi,x = GGI$unrate1,p_pseudo = 2)
selection.lag.lead(GGI$stdcpi,x = GGI$unrate2,p_pseudo = 10)

GGI$MAR$table2 = matrix(0,3,2)
GGI$MAR$table2[,1] = GGI$MAR$select$loglikelihood
GGI$MAR$table2[,2] = c(Box.test(GGI$MAR$model11$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGI$MAR$model12$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGI$MAR$model13$residuals,lag = 12,type = "Ljung-Box")$p.value)

###Great Moderation Period###

n = length(GGM$cpi)
nlag = 12
GGM$lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        GGM$cpi[1:(n - x)]))
#lags = lags[-c(1:nlag),]

GGM$AIC = c()

for(i in 1:nlag){
  exo = GGM$unrate1[-c(1:i)]
  auxm = lm(GGM$cpi[-c(1:i)]~1+GGM$lags[-c(1:i),1:i]+exo)
  GGM$AIC[i] = AIC(auxm)}


GGM$MAR$model11 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 2,p_NC = 0)
GGM$MAR$model12 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 1,p_NC = 1)
GGM$MAR$model13 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 0,p_NC = 2)

GGM$MAR$table = matrix(0,6,8)
GGM$MAR$table[1,]= c(GGM$MAR$model11$coefficients[c(1,2,3,4,4,5,6,7)])
GGM$MAR$table[2,]= c(GGM$MAR$model11$se[c(1,2,3,4,4,5,6,7)])
GGM$MAR$table[3,]= c(GGM$MAR$model12$coefficients[c(1,2,2,3,3,4,5,6)])
GGM$MAR$table[4,]= c(GGM$MAR$model12$se[c(1,2,2,3,3,4,5,6)])
GGM$MAR$table[5,]= c(GGM$MAR$model13$coefficients[c(1,2,2,3,4,5,6,7)])
GGM$MAR$table[6,]= c(GGM$MAR$model13$se[c(1,2,2,3,4,5,6,7)])
rownames(GGM$MAR$table) = c("AR(2,0)","s.e.", "MAR(1,1)","s.e.", "AR(0,2)","s.e.")
colnames(GGM$MAR$table) = c("int", "lag1", "lag2", "lead1", "lead2", "exo", "df", "scale" )

GGM$MAR$select = selection.lag.lead(GGM$cpi,x = GGM$unrate1,p_pseudo = 2)
selection.lag.lead(GGM$cpi,x = GGM$unrate1,p_pseudo = 12)

GGM$MAR$select5 = selection.lag.lead(GGM$cpi,x = GGM$unrate1,p_pseudo = 5)

GGM$MAR$model50 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 5,p_NC = 0)
GGM$MAR$model41 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 4,p_NC = 1)
GGM$MAR$model32 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 3,p_NC = 2)
GGM$MAR$model23 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 2,p_NC = 3)
GGM$MAR$model14 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 1,p_NC = 4)
GGM$MAR$model05 = mixed(y = GGM$cpi,x = GGM$unrate1,p_C = 0,p_NC = 5)

c(Box.test(GGM$MAR$model50$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GGM$MAR$model41$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GGM$MAR$model32$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GGM$MAR$model23$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GGM$MAR$model14$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GGM$MAR$model05$residuals,lag = 12,type = "Ljung-Box")$p.value)


GGM$MAR$table2 = matrix(0,3,2)
GGM$MAR$table2[,1] = GGM$MAR$select2$loglikelihood
GGM$MAR$table2[,2] = c(Box.test(GGM$MAR$model11$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGM$MAR$model12$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGM$MAR$model13$residuals,lag = 12,type = "Ljung-Box")$p.value)

GGM$MAR$table3 = matrix(0,6,2)
GGM$MAR$table3[,1] = GGM$MAR$select5$loglikelihood
GGM$MAR$table3[,2] = c(Box.test(GGM$MAR$model50$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGM$MAR$model41$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGM$MAR$model32$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGM$MAR$model23$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGM$MAR$model14$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGM$MAR$model05$residuals,lag = 12,type = "Ljung-Box")$p.value)


###Great Recession Period###

n = length(GGR$cpi)
nlag = 12
GGR$lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        GGR$cpi[1:(n - x)]))
#lags = lags[-c(1:nlag),]

GGR$AIC = c()

for(i in 1:nlag){
  exo = GGR$unrate1[-c(1:i)]
  auxm = lm(GGR$cpi[-c(1:i)]~1+GGR$lags[-c(1:i),1:i]+exo)
  GGR$AIC[i] = AIC(auxm)}

mixed(y = GGR$cpi,x = GGR$unrate1,p_C = 0,p_NC = 1)
mixed(y = GGR$cpi,x = GGR$unrate1,p_C = 1,p_NC = 0)

GGR$MAR$model11 = mixed(y = GGR$cpi,x = GGR$unrate1,p_C = 2,p_NC = 0)
GGR$MAR$model12 = mixed(y = GGR$cpi,x = GGR$unrate1,p_C = 1,p_NC = 1)
GGR$MAR$model13 = mixed(y = GGR$cpi,x = GGR$unrate1,p_C = 0,p_NC = 2)

GGR$MAR$table = matrix(0,6,8)
GGR$MAR$table[1,]= c(GGR$MAR$model11$coefficients[c(1,2,3,4,4,5,6,7)])
GGR$MAR$table[2,]= c(GGR$MAR$model11$se[c(1,2,3,4,4,5,6,7)])
GGR$MAR$table[3,]= c(GGR$MAR$model12$coefficients[c(1,2,2,3,3,4,5,6)])
GGR$MAR$table[4,]= c(GGR$MAR$model12$se[c(1,2,2,3,3,4,5,6)])
GGR$MAR$table[5,]= c(GGR$MAR$model13$coefficients[c(1,2,2,3,4,5,6,7)])
GGR$MAR$table[6,]= c(GGR$MAR$model13$se[c(1,2,2,3,4,5,6,7)])
rownames(GGR$MAR$table) = c("AR(2,0)","s.e.", "MAR(1,1)","s.e.", "AR(0,2)","s.e.")
colnames(GGR$MAR$table) = c("int", "lag1", "lag2", "lead1", "lead2", "exo", "df", "scale" )


GGR$MAR$select = selection.lag.lead(GGR$cpi,x = GGR$unrate1,p_pseudo = 2)

GGR$MAR$table2 = matrix(0,3,2)
GGR$MAR$table2[,1] = GGR$MAR$select$loglikelihood
GGR$MAR$table2[,2] = c(Box.test(GGR$MAR$model11$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGR$MAR$model12$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GGR$MAR$model13$residuals,lag = 12,type = "Ljung-Box")$p.value)







p = c(0,0,0,0,0)
#check = optim(p,LAD,yt=GGI$cpi,r=NULL,order = c(1,1),hessian=TRUE)
check2 = mixed(y = GGI$stdcpi, x = NULL, p_C = 1, p_NC = 1)
p = c(c(check$par),1,5)
check3 = optim(p,TMLE,yt=GGI$stdcpi,r=NULL,order = c(1,1),hessian=TRUE,
               control = list(maxit = 5000))
nlag = 2
nlead = 2
n = length(GGI$stdcpi)
auxLA =sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                    GGI$stdcpi[1:(n - x)]))
auxLE = sapply(1:nlead, function(x) c(GGI$stdcpi[-c(1:x)],
                                      rep(NA, length.out = x)))

aux = cbind(auxLE,auxLA)

p = rep(0,4)
nstd = length(GGI$stdcpi)
p = rep(0,4)
check4 = optim(p,LAD22,yt=GGI$stdcpi[-c(1,2,nstd-1,nstd)],
               r=aux[-c(1,2,nstd-1,nstd),],
               hessian=TRUE,
               control = list(maxit = 5000))

LAD22(p = check4$par,yt = GGI$stdcpi[-c(1,2,nstd-1,nstd)],r=aux[-c(1,2,nstd-1,nstd),])
LAD2(p = c(0,check4$par),yt = GGI$stdcpi,r = NULL,order = c(2,2))

p = rep(0,6)
check5 = optim(par = p,fn = LAD2,yt = GGI$stdcpi, r = GGI$unrate1,order = c(2,2),
               control = list(maxit = 2000),hessian = TRUE)
mixed(y = GGI$stdcpi,x = NULL,p_C = 2,p_NC = 2)$coefficients

optim(c(rep(0,4),4,4),t22,yt=GGI$stdcpi[-c(1,2,nstd-1,nstd)],
      r=aux[-c(1,2,nstd-1,nstd),],
      hessian=TRUE,
      control = list(maxit = 5000))$par
pp = rep(0.5,5)
LAD2(p = c(0,pp),yt = GGI$stdcpi,r = GGI$unrate1,order = c(1,1))
LAD2(p = c(0,pp),yt = GGI$stdcpi,order = c(2,2))

LAD11(p = pp,yt = GGI$stdcpi[-c(1,nstd)],r = aux[-c(1,nstd),c(3,1)])
LAD22(p = pp,yt = GGI$stdcpi[-c(1,2,nstd-1,nstd)],r = aux[-c(1,2,nstd-1,nstd),])
LADest(yt = GGI$stdcpi,exoreg = GGI$unrate1,order = c(3,1))$coefficients
t(mixed(y = GGI$stdcpi,x = GGI$unrate1,p_C = 3,p_NC = 1)$coefficients)

TMLE(p = c(0,pp,3,3), yt= GGI$stdcpi,r = GGI$unrate1,order = c(1,1))
TMLE(p = c(0,pp,3,3), yt= GGI$stdcpi,r = GGI$unrate1,order = c(3,1))

TMLEest(yt = GGI$stdcpi,exoreg = rep(0,length(GGI$stdcpi)),order = c(3,1))$coefficients



LAD22=function(p, yt, r){
  const = p[1]
  p = p[-1]
  cond=yt-p[1]*r[,1]-p[2]*r[,2]-p[3]*r[,3]+
    p[1]*p[3]*yt+p[2]*p[3]*r[,1]-
    p[4]*r[,4]+p[1]*p[4]*r[,3]+p[2]*p[4]*yt
  cond = cond - const
  stat=sum(abs(cond))
  return(stat)
}

t22e=function(p,yt,r, exo){
  const = p[1]
  p = p[-1]
  sigma=p[5]
  lambda=p[6]
  pexo = p[-c(1:6)]
  cond=yt-p[1]*r[,1]-p[2]*r[,2]-p[3]*r[,3]+
    p[1]*p[3]*yt+p[2]*p[3]*r[,1]-
    p[4]*r[,4]+p[1]*p[4]*r[,3]+p[2]*p[4]*yt
  cond = cond - const - pexo%*%t(exo)
  ll1=log(gamma((lambda+1)/2)/gamma(lambda/2)/sqrt(pi))-0.5*log(lambda-2)-log(sigma)
  ll2=-(lambda+1)/2*log(1+cond^2/sigma^2/(lambda-2))
  loglike=-(sum(ll1+ll2))
}

pt = c(pp,3,3)
LAD22(p = pp,yt = GGI$stdcpi[-c(1,2,nstd-1,nstd)],r = aux[-c(1,2,nstd-1,nstd),])
optim(par = pp,LAD22,yt = GGI$stdcpi[-c(1,2,nstd-1,nstd)],
      r = aux[-c(1,2,nstd-1,nstd),], control = list(maxit = 1000))
optim(par = pt,t22,yt = GGI$stdcpi[-c(1,2,nstd-1,nstd)],
      r = aux[-c(1,2,nstd-1,nstd),], control = list(maxit = 1000))
mixed(y = GGI$stdcpi,x = NULL,p_C = 2,p_NC = 2)$coefficients

pte = c(pt,0)
optim(par = pte,t22e,yt = GGI$stdcpi[-c(1,2,nstd-1,nstd)],
      r = aux[-c(1,2,nstd-1,nstd),], 
      exo = GGI$unrate1[-c(1,2,nstd-1,nstd)],control = list(maxit = 2000))
mixed(y = GGI$stdcpi,x = GGI$unrate1,p_C = 2,p_NC = 2)$coefficients



