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

#hpfiltered = hpfilter(x = unrate,type = "lambda",freq = 14400)
#hpfiltered2 = hpfilter(x = unrate,type = "lambda",freq = 28800)


#ts.plot(unrate)
#lines(hpfiltered$trend, col = "red")
#lines(hpfiltered$cycle+6, col = "blue")
#abline(a=6,b = 0)
#var(hpfiltered$cycle)

# ts.plot(unrate)
# lines(hpfiltered2$trend, col = "red")
# lines(hpfiltered2$cycle+6, col = "blue")
# abline(a=6,b = 0)
# var(hpfiltered2$cycle)
#adf.test(hpfiltered$cycle)

###Baxter-King filter####

#bkfiltered1 = bkfilter(unrate,pl = 2,pu = 96,nfix = 36) 
#filter out the irregular component and the business cycle component

bkfiltered2 = bkfilter(unrate,pl = 18,pu = 96,nfix = 36) 
#filter out the business cycle of 1.5-8 years

#bkfiltered11 = bkfilter(unrate,pl = 2,pu = 96,nfix = 24) 
#filter out the business cycle of 1.5-8 years

bkfiltered12 = bkfilter(unrate,pl = 18,pu = 96,nfix = 24) 
#filter out the business cycle of 1.5-8 years


#ts.plot(bkfiltered1$cycle, col = "red")
#lines(bkfiltered12$cycle, col = "blue")
#lines(bkfiltered$cycle, col = "green4")
#abline(0,0)

#ts.plot(bkfiltered2$cycle, col = "red")
#lines(bkfiltered12$cycle, col = "blue")
#lines(bkfiltered$cycle, col = "green4")
#abline(0,0)

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

corecpi = data.frame(read_excel("CPILFESL_M1957-2019.xls",sheet = "data"))[,2]
corecpi_log = ts(data = diff(log(corecpi)),start = c(1957,1),end = c(2019,5),frequency = 12)
corecpi.ts = window(x = corecpi_log,start = c(1957,1),end = c(2019,4))

ppi = data.frame(read_excel("PPIACO_M1913-2019.xls",sheet = "data"))[,2]
ppi_log = ts(data = diff(log(ppi)),start = c(1913,1),end = c(2019,5),frequency = 12)
ppi.ts = window(x = ppi_log,start = c(1953,1),end = c(2019,4))

unratecycle.ts = ts(data = bkfiltered12$cycle,start = c(1953,1),end = c(2019,5),frequency = 12)
#unratecycle2.ts = ts(data = bkfiltered11$cycle,start = c(1953,1),end = c(2019,5),frequency = 12)

##Above we split the data in appropriate time range matching the unemploymen sequene if possible

cpi1.ts = cpi.ts[!is.na(unratecycle.ts[-length(unratecycle.ts)])]
#corecpi1.ts = corecpi.ts[!is.na(unratecycle.ts[-length(unratecycle.ts)])]
ppi1.ts = ppi.ts[!is.na(unratecycle.ts[-length(unratecycle.ts)])]

#cpi1.gr = cpi.grate[!is.na(unratecycle.ts[-length(unratecycle.ts)])]

unratecycle.ts = unratecycle.ts[!is.na(unratecycle.ts)]
#unratecycle2.ts = unratecycle2.ts[!is.na(unratecycle2.ts)]

#ts.plot(cpi.ts*100)
#lines(unratecycle.ts, col = "red")
#lines(unratecycle2.ts, col = "blue")


unratecycle.ts = ts(data = unratecycle.ts,start = c(1955,1),end = c(2017,5),frequency = 12)
#unratecycle2.ts = ts(data = unratecycle2.ts,start = c(1955,1),end = c(2017,5),frequency = 12)

cpi2.ts = ts(data = cpi1.ts,start = c(1955,1),end = c(2017,5),frequency = 12)
corecpi2.ts = window(x = corecpi.ts,start = c(1957,1),end = c(2017,5),frequency = 12)
ppi2.ts = ts(data = ppi1.ts,start = c(1955,1),end = c(2017,5),frequency = 12)

#cpi2.gr = ts(data = cpi1.gr,start = c(1955,1),end = c(2017,5),frequency = 12)


###Divide the data into three major periods
##Great Inflation: 1960 - 1983
GI = list()
GI$cpi = window(x = cpi2.ts, start = c(1965,1), end = c(1982,12))
GI$corecpi= window(x = corecpi2.ts, start = c(1965,1), end = c(1982,12))
GI$ppi = window(x = ppi2.ts, start = c(1965,1), end = c(1982,12))

#GI$cpi = GI$cpi2

#GI$unrate1 = window(x = unratecycle.ts, start = c(1960,1), end = c(1982,12))
GI$unrate = window(x = unratecycle.ts, start = c(1965,1), end = c(1982,12))

#GI$unrate = GI$unrate3
GI$n = length(GI$cpi)
#GI$stdcpi = (GI$cpi - mean(GI$cpi))/sqrt(var(GI$cpi))


#check = modelfit2(GI$cpi,l = 20)
GI$lagselect = modelfit2(GI$cpi,l = 20,exo = GI$unrate)
GI$core$lagselect = modelfit2(GI$corecpi,l = 20,exo = GI$unrate)
GI$PPIN$lagselect = modelfit2(GI$ppi,l = 20,exo = GI$unrate)


GI$unitroot0 = summary(ur.df(GI$cpi,type = "none"))
GI$unitrootD = summary(ur.df(GI$cpi,type = "drift"))
GI$unitrootT = summary(ur.df(GI$cpi,type = "trend"))

GI$core$unitroot0 = summary(ur.df(GI$corecpi,type = "none"))
GI$core$unitrootD = summary(ur.df(GI$corecpi,type = "drift"))
GI$core$unitrootT = summary(ur.df(GI$corecpi,type = "trend"))

GI$PPIN$unitroot0 = summary(ur.df(GI$ppi,type = "none"))
GI$PPIN$unitrootD = summary(ur.df(GI$ppi,type = "drift"))
GI$PPIN$unitrootT = summary(ur.df(GI$ppi,type = "trend"))


GI$causalmodel = arima(x = GI$cpi,order = c(2,0,0),
                        xreg = GI$unrate,method = "CSS")
GI$core$causalmodel = arima(x = GI$corecpi,order = c(2,0,0),
                        xreg = GI$unrate,method = "CSS")
GI$PPIN$causalmodel = arima(x = GI$ppi,order = c(1,0,0),
                             xreg = GI$unrate,method = "CSS")

GI$ACtest1 = diagnostics(cpi = GI$cpi, unrate = GI$unrate,range = 20,order = 30)
GI$ACtest2 = diagnostics(cpi = GI$cpi, unrate = GI$unrate,range = 20,order = 10)

GI$core$ACtest1 = diagnostics(cpi = GI$corecpi, unrate = GI$unrate,range = 20,order = 30)
GI$core$ACtest2 = diagnostics(cpi = GI$corecpi, unrate = GI$unrate,range = 20,order = 10)

GI$PPIN$ACtest1 = diagnostics(cpi = GI$ppi, unrate = GI$unrate,range = 20,order = 30)
GI$PPIN$ACtest2 = diagnostics(cpi = GI$ppi, unrate = GI$unrate,range = 20,order = 10)
GI$PPIN$ACtest3 = diagnostics(cpi = GI$ppi, unrate = GI$unrate,range = 20,order = 20)


ts.plot(GI$ACtest1$BGp.value)
which(GI$ACtest2$BGp.value<0.05)#9 lags for both orders

ts.plot(GI$core$ACtest1$BGp.value)
which(GI$core$ACtest2$BGp.value<0.05) #6 lags for order 10 or 9 lags for order 30

ts.plot(GI$PPIN$ACtest3$BGp.value)
which(GI$PPIN$ACtest3$BGp.value<0.05) #11 lags for order 10 and 3 lags for order 30

##Great Moderation: 1985 - 2007
GM = list()
GM$cpi = window(x = cpi2.ts, start = c(1983,1), end = c(2006,12))
GM$corecpi = window(x = corecpi2.ts, start = c(1983,1), end = c(2006,12))
GM$ppi = window(x = ppi2.ts, start = c(1983,1), end = c(2006,12))

GM$unrate = window(x = unratecycle.ts, start = c(1983,1), end = c(2006,12))
#GM$unrate2 = window(x = unratecycle2.ts, start = c(1983,1), end = c(2006,12))
GM$n = length(GM$cpi)
#GM$stdcpi = (GM$cpi - mean(GM$cpi))/sqrt(var(GM$cpi))

GM$lagselect = modelfit2(GM$cpi,l = 20,exo = GM$unrate)
GM$core$lagselect = modelfit2(GM$corecpi,l = 20,exo = GM$unrate)
GM$PPIN$lagselect = modelfit2(GM$ppi,l = 20,exo = GM$unrate)

#GM$lagselect2 = modelfit2(GM$cpi,l = 20,exo = GM$unrate2)

GM$ACtest1 = diagnostics(cpi = GM$cpi, unrate = GM$unrate,range = 20,order = 30)
GM$ACtest2 = diagnostics(cpi = GM$cpi, unrate = GM$unrate,range = 20,order = 10)

GM$core$ACtest1 = diagnostics(cpi = GM$corecpi, unrate = GM$unrate,range = 20,order = 30)
GM$core$ACtest2 = diagnostics(cpi = GM$corecpi, unrate = GM$unrate,range = 20,order = 10)

GM$PPIN$ACtest1 = diagnostics(cpi = GM$ppi, unrate = GM$unrate,range = 20,order = 30)
GM$PPIN$ACtest2 = diagnostics(cpi = GM$ppi, unrate = GM$unrate,range = 20,order = 10)

GM$unitroot0 = summary(ur.df(GM$cpi,type = "none"))
GM$unitrootD = summary(ur.df(GM$cpi,type = "drift"))
GM$unitrootT = summary(ur.df(GM$cpi,type = "trend"))

GM$core$unitroot0 = summary(ur.df(GM$corecpi,type = "none"))
GM$core$unitrootD = summary(ur.df(GM$corecpi,type = "drift"))
GM$core$unitrootT = summary(ur.df(GM$corecpi,type = "trend"))

GM$PPIN$unitroot0 = summary(ur.df(GM$ppi,type = "none"))
GM$PPIN$unitrootD = summary(ur.df(GM$ppi,type = "drift"))
GM$PPIN$unitrootT = summary(ur.df(GM$ppi,type = "trend"))

GM$causalmodel = arima(x = GM$cpi,order = c(2,0,0),
                        xreg = GM$unrate,method = "CSS")

GM$core$causalmodel = arima(x = GM$corecpi,order = c(6,0,0),
                        xreg = GM$unrate,method = "CSS")

GM$PPIN$causalmodel = arima(x = GM$ppi,order = c(1,0,0),
                        xreg = GM$unrate,method = "CSS")

##Great Recession: 2007 - 2017
GR = list()
GR$cpi = window(x = cpi2.ts, start = c(2007,1), end = c(2017,5))
GR$corecpi = window(x = corecpi2.ts, start = c(2007,1), end = c(2017,5))
GR$ppi = window(x = ppi2.ts, start = c(2007,1), end = c(2017,5))

GR$unrate = window(x = unratecycle.ts, start = c(2007,1), end = c(2017,5))
#GR$unrate2 = window(x = unratecycle2.ts, start = c(2007,1), end = c(2017,5))

GR$n = length(GR$cpi)

#GR$stdcpi = (GR$cpi - mean(GR$cpi))/sqrt(var(GR$cpi))

GR$lagselect = modelfit2(GR$cpi,l = 20,exo = GR$unrate)
GR$core$lagselect = modelfit2(GR$corecpi,l = 20,exo = GR$unrate)
GR$PPIN$lagselect = modelfit2(GR$ppi,l = 20,exo = GR$unrate)


#GR$lagselect2 = modelfit2(GR$cpi,l = 20,exo = GR$unrate2)

GR$ACtest1 = diagnostics(cpi = GR$cpi, unrate = GR$unrate,range = 20,order = 30)
GR$ACtest2 = diagnostics(cpi = GR$cpi, unrate = GR$unrate,range = 20,order = 10)

GR$core$ACtest1 = diagnostics(cpi = GR$corecpi, unrate = GR$unrate,range = 20,order = 30)
GR$core$ACtest2 = diagnostics(cpi = GR$corecpi, unrate = GR$unrate,range = 20,order = 10)

GR$PPIN$ACtest1 = diagnostics(cpi = GR$ppi, unrate = GR$unrate,range = 20,order = 30)
GR$PPIN$ACtest2 = diagnostics(cpi = GR$ppi, unrate = GR$unrate,range = 20,order = 10)

GR$unitroot0 = summary(ur.df(GR$cpi,type = "none"))
GR$unitrootD = summary(ur.df(GR$cpi,type = "drift"))
GR$unitrootT = summary(ur.df(GR$cpi,type = "trend"))

GR$core$unitroot0 = summary(ur.df(GR$corecpi,type = "none"))
GR$core$unitrootD = summary(ur.df(GR$corecpi,type = "drift"))
GR$core$unitrootT = summary(ur.df(GR$corecpi,type = "trend"))

GR$PPIN$unitroot0 = summary(ur.df(GR$ppi,type = "none"))
GR$PPIN$unitrootD = summary(ur.df(GR$ppi,type = "drift"))
GR$PPIN$unitrootT = summary(ur.df(GR$ppi,type = "trend"))

GR$causalmodel = arima(x = GR$cpi,order = c(2,0,0),
                        xreg = GR$unrate,method = "CSS")

GR$core$causalmodel = arima(x = GR$corecpi,order = c(2,0,0),
                        xreg = GR$unrate,method = "CSS")

GR$PPIN$causalmodel = arima(x = GR$ppi,order = c(3,0,0),
                        xreg = GR$unrate,method = "CSS")

####Tables####
library(e1071)  
tables = list()
#Summary statistics table CPI
tables$sumstat = matrix(0,10,3)
tables$sumstat[,1] = c(mean(GI$cpi),median(GI$cpi),var(GI$cpi),skewness(GI$cpi),
                       kurtosis(GI$cpi),min(GI$cpi),max(GI$cpi),
                       GI$unitroot0@teststat,GI$unitrootD@teststat[1],
                       GI$unitrootT@teststat[1])
tables$sumstat[,2] = c(mean(GM$cpi),median(GM$cpi),var(GM$cpi),skewness(GM$cpi),
                       kurtosis(GM$cpi),min(GM$cpi),max(GM$cpi),
                       GM$unitroot0@teststat,GM$unitrootD@teststat[1],
                       GM$unitrootT@teststat[1])
tables$sumstat[,3] = c(mean(GR$cpi),median(GR$cpi),var(GR$cpi),skewness(GR$cpi),
                       kurtosis(GR$cpi),min(GR$cpi),max(GR$cpi),
                       GR$unitroot0@teststat,GR$unitrootD@teststat[1],
                       GR$unitrootT@teststat[1])
rownames(tables$sumstat) = c("Mean","Median","Variance","Skewness","Kurtosis",
                             "Min", "Max", "ADF test (origin)", 
                             "ADF test (const)","ADF test (const+trend)")

#Summary statistics table CPI core
tables$core$sumstat = matrix(0,10,3)
tables$core$sumstat[,1] = c(mean(GI$corecpi),median(GI$corecpi),var(GI$corecpi),skewness(GI$corecpi),
                       kurtosis(GI$corecpi),min(GI$corecpi),max(GI$corecpi),
                       GI$core$unitroot0@teststat,GI$core$unitrootD@teststat[1],
                       GI$core$unitrootT@teststat[1])
tables$core$sumstat[,2] = c(mean(GM$corecpi),median(GM$corecpi),var(GM$corecpi),skewness(GM$corecpi),
                       kurtosis(GM$corecpi),min(GM$corecpi),max(GM$corecpi),
                       GM$core$unitroot0@teststat,GM$core$unitrootD@teststat[1],
                       GM$core$unitrootT@teststat[1])
tables$core$sumstat[,3] = c(mean(GR$corecpi),median(GR$corecpi),var(GR$corecpi),skewness(GR$corecpi),
                       kurtosis(GR$corecpi),min(GR$corecpi),max(GR$corecpi),
                       GR$core$unitroot0@teststat,GR$core$unitrootD@teststat[1],
                       GR$core$unitrootT@teststat[1])
rownames(tables$core$sumstat) = c("Mean","Median","Variance","Skewness","Kurtosis",
                             "Min", "Max", "ADF test (origin)", 
                             "ADF test (const)","ADF test (const+trend)")

#Summary statistics table PPI
tables$PPIN$sumstat = matrix(0,10,3)
tables$PPIN$sumstat[,1] = c(mean(GI$ppi),median(GI$ppi),var(GI$ppi),skewness(GI$ppi),
                            kurtosis(GI$ppi),min(GI$ppi),max(GI$ppi),
                            GI$PPIN$unitroot0@teststat,GI$PPIN$unitrootD@teststat[1],
                            GI$PPIN$unitrootT@teststat[1])
tables$PPIN$sumstat[,2] = c(mean(GM$ppi),median(GM$ppi),var(GM$ppi),skewness(GM$ppi),
                            kurtosis(GM$ppi),min(GM$ppi),max(GM$ppi),
                            GM$PPIN$unitroot0@teststat,GM$PPIN$unitrootD@teststat[1],
                            GM$PPIN$unitrootT@teststat[1])
tables$PPIN$sumstat[,3] = c(mean(GR$ppi),median(GR$ppi),var(GR$ppi),skewness(GR$ppi),
                            kurtosis(GR$ppi),min(GR$ppi),max(GR$ppi),
                            GR$PPIN$unitroot0@teststat,GR$PPIN$unitrootD@teststat[1],
                            GR$PPIN$unitrootT@teststat[1])
rownames(tables$PPIN$sumstat) = c("Mean","Median","Variance","Skewness","Kurtosis",
                                  "Min", "Max", "ADF test (origin)", 
                                  "ADF test (const)","ADF test (const+trend)")



#Diagnostics test table
tables$tests = matrix(0,18,21)
tables$tests[1:6,] = rbind(GI$ACtest1$LBtest,GI$ACtest1$LBp.value,
                           GI$ACtest1$BGtest,GI$ACtest1$BGp.value,
                           GI$ACtest1$JBstat,GI$ACtest1$JBp.value)
tables$tests[7:12,] = rbind(GM$ACtest1$LBtest,GM$ACtest1$LBp.value,
                            GM$ACtest1$BGtest,GM$ACtest1$BGp.value,
                            GM$ACtest1$JBstat,GM$ACtest1$JBp.value)
tables$tests[13:18,] = rbind(GR$ACtest1$LBtest,GR$ACtest1$LBp.value,
                             GR$ACtest1$BGtest,GR$ACtest1$BGp.value,
                             GR$ACtest1$JBstat,GR$ACtest1$JBp.value)
tables$testrownames = c("Ljung-Box, ord=30", "L-B p.value",
                        "Breusch-Godfrey, ord=30", "B-G p.value",
                        "Jarque-Bera", "J-B p.value")
rownames(tables$tests) = rep(tables$testrownames,3)
colnames(tables$tests) = c(1:21)

#diagnostics tests table for order 10
tables$tests2 = matrix(0,18,21)
tables$tests2[1:6,] = rbind(GI$ACtest2$LBtest,GI$ACtest2$LBp.value,
                            GI$ACtest2$BGtest,GI$ACtest2$BGp.value,
                            GI$ACtest2$JBstat,GI$ACtest2$JBp.value)
tables$tests2[7:12,] = rbind(GM$ACtest2$LBtest,GM$ACtest2$LBp.value,
                             GM$ACtest2$BGtest,GM$ACtest2$BGp.value,
                             GM$ACtest2$JBstat,GM$ACtest2$JBp.value)
tables$tests2[13:18,] = rbind(GR$ACtest2$LBtest,GR$ACtest2$LBp.value,
                              GR$ACtest2$BGtest,GR$ACtest2$BGp.value,
                              GR$ACtest2$JBstat,GR$ACtest2$JBp.value)
tables$testrownames = c("Ljung-Box, ord=10", "L-B p.value",
                        "Breusch-Godfrey, ord=10", "B-G p.value",
                        "Jarque-Bera", "J-B p.value")
rownames(tables$tests2) = rep(tables$testrownames,3)
colnames(tables$tests2) = c(1:21)

#AICBIC table
tables$aicbic = matrix(0,6,20)
tables$aicbic[1:2,] = rbind(GI$lagselect$aic,GI$lagselect$bic)
tables$aicbic[3:4,] = rbind(GM$lagselect$aic,GM$lagselect$bic)
tables$aicbic[5:6,] = rbind(GR$lagselect$aic,GR$lagselect$bic)
rownames(tables$aicbic) = rep(c("AIC","BIC"),3)


####Step 3. Analysis of the noncausal setings####
###Follow Lanne, Saikonen lag order selection

###Great Inflation Period###
n = length(GI$cpi)
nlag = 12
GI$lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        GI$cpi[1:(n - x)]))
#lags = lags[-c(1:nlag),]

GI$AIC = c()
GI$AIC2 = c()

for(i in 1:nlag){
  exo = GI$unrate[-c(1:i)]
  auxm = lm(GI$cpi[-c(1:i)]~1+GI$lags[-c(1:i),1:i]+exo)
  GI$AIC[i] = AIC(auxm)}

# GI$lags2 = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
#                                          GI$stdcpi[1:(n - x)]))
# for(i in 1:nlag){
#   exo = GI$unrate2[-c(1:i)]
#   auxm = lm(GI$stdcpi[-c(1:i)]~1+GI$lags2[-c(1:i),1:i]+exo)
#   GI$AIC2[i] = AIC(auxm)}

#Lag selection process is a sequential procedure
#therefore simple commands from MARX package are not appropriate

#=>three possible candidates c(2,5,9)

GI$lagorder = c(2,5,9)
exo2 = GI$unrate[-c(1:GI$lagorder[1])]
GI$cm2 = lm(GI$cpi[-c(1:GI$lagorder[1])]~1+
               GI$lags[-c(1:GI$lagorder[1]),1:GI$lagorder[1]]+exo2)


exo5 = GI$unrate[-c(1:GI$lagorder[2])]
GI$cm5 = lm(GI$cpi[-c(1:GI$lagorder[2])]~1+
               GI$lags[-c(1:GI$lagorder[2]),1:GI$lagorder[2]]+exo5)

exo9 = GI$unrate[-c(1:GI$lagorder[3])]
GI$cm9 = lm(GI$cpi[-c(1:GI$lagorder[3])]~1+
               GI$lags[-c(1:GI$lagorder[3]),1:GI$lagorder[3]]+exo9)


GI$cmstat2 = summary(GI$cm2)
GI$cmstat5 = summary(GI$cm5)
GI$cmstat9 = summary(GI$cm9)

GI$MAR = list()
GI$MAR$causal2 = selection.lag.lead(y = GI$cpi, x = GI$unrate,p_pseudo = 2)

GI$MAR$causal5 = selection.lag.lead(y = GI$cpi,x = GI$unrate,p_pseudo = 5)

GI$MAR$causal9 = selection.lag.lead(y = GI$cpi,x = GI$unrate,p_pseudo = 9)

GI$MAR$model11 = mixed(y = GI$cpi,x = GI$unrate,p_C = 2,p_NC = 0)
GI$MAR$model12 = mixed(y = GI$cpi,x = GI$unrate,p_C = 1,p_NC = 1)
GI$MAR$model13 = mixed(y = GI$cpi,x = GI$unrate,p_C = 0,p_NC = 2)

GI$core$MAR$causal2 = selection.lag.lead(y = GI$corecpi, x = GI$unrate,p_pseudo = 6)


#GM$core$MAR$causal2 = selection.lag.lead(y = GM$corecpi, x = GM$unrate,p_pseudo = 6)
#GR$core$MAR$causal2 = selection.lag.lead(y = GR$corecpi, x = GR$unrate,p_pseudo = 2)



GI$core$MAR$model11 = mixed(y = GI$corecpi,x = GI$unrate,p_C = 2,p_NC = 0)
GI$core$MAR$model12 = mixed(y = GI$corecpi,x = GI$unrate,p_C = 1,p_NC = 1)
GI$core$MAR$model13 = mixed(y = GI$corecpi,x = GI$unrate,p_C = 0,p_NC = 2)


GI$MAR$table = matrix(0,6,8)
GI$MAR$table[1,]= c(GI$MAR$model11$coefficients[c(1,2,3,4,4,5,6,7)])
GI$MAR$table[2,]= c(GI$MAR$model11$se[c(1,2,3,4,4,5,6,7)])
GI$MAR$table[3,]= c(GI$MAR$model12$coefficients[c(1,2,2,3,3,4,5,6)])
GI$MAR$table[4,]= c(GI$MAR$model12$se[c(1,2,2,3,3,4,5,6)])
GI$MAR$table[5,]= c(GI$MAR$model13$coefficients[c(1,2,2,3,4,5,6,7)])
GI$MAR$table[6,]= c(GI$MAR$model13$se[c(1,2,2,3,4,5,6,7)])
rownames(GI$MAR$table) = c("AR(2,0)","s.e.", "MAR(1,1)","s.e.", "AR(0,2)","s.e.")
colnames(GI$MAR$table) = c("int", "lag1", "lag2", "lead1", "lead2", "exo", "df", "scale" )

GI$MAR$model11$MSE = sum((GI$cpi[-c(1,2)]-GI$cm2$fitted.values)^2)
GI$MAR$model12$MSE = sum((GI$cpi[-c(1,n)]-GI$MAR$model12$fitted.values)^2)
GI$MAR$model13$MSE = sum((GI$cpi[-c(n-1,n)]-GI$MAR$model13$fitted.values)^2)

nlead = 2
GI$leads = sapply(1:nlead, function(x) c(GI$cpi[-c(1:x)],rep(NA, length.out = x)))


Box.test(GI$MAR$model11$residuals,lag = 12,type = "Ljung-Box")
Box.test(GI$MAR$model12$residuals,lag = 12,type = "Ljung-Box")
Box.test(GI$MAR$model13$residuals,lag = 12,type = "Ljung-Box")

selection.lag.lead(GI$cpi,x = GI$unrate,p_pseudo = 2)
GI$MAR$select = selection.lag.lead(GI$cpi,x = GI$unrate,p_pseudo = 2)

# selection.lag.lead(GI$stdcpi,x = GI$unrate1,p_pseudo = 2)
# selection.lag.lead(GI$stdcpi,x = GI$unrate2,p_pseudo = 10)

GI$MAR$table2 = matrix(0,3,2)
GI$MAR$table2[,1] = GI$MAR$select$loglikelihood
GI$MAR$table2[,2] = c(Box.test(GI$MAR$model11$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GI$MAR$model12$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GI$MAR$model13$residuals,lag = 12,type = "Ljung-Box")$p.value)

###Great Moderation Period###

n = length(GM$cpi)
nlag = 12
GM$lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        GM$cpi[1:(n - x)]))
#lags = lags[-c(1:nlag),]

GM$AIC = c()

for(i in 1:nlag){
  exo = GM$unrate[-c(1:i)]
  auxm = lm(GM$cpi[-c(1:i)]~1+GM$lags[-c(1:i),1:i]+exo)
  GM$AIC[i] = AIC(auxm)}


GM$MAR$model11 = mixed(y = GM$cpi,x = GM$unrate,p_C = 2,p_NC = 0)
GM$MAR$model12 = mixed(y = GM$cpi,x = GM$unrate,p_C = 1,p_NC = 1)
GM$MAR$model13 = mixed(y = GM$cpi,x = GM$unrate,p_C = 0,p_NC = 2)

GM$MAR$table = matrix(0,6,8)
GM$MAR$table[1,]= c(GM$MAR$model11$coefficients[c(1,2,3,4,4,5,6,7)])
GM$MAR$table[2,]= c(GM$MAR$model11$se[c(1,2,3,4,4,5,6,7)])
GM$MAR$table[3,]= c(GM$MAR$model12$coefficients[c(1,2,2,3,3,4,5,6)])
GM$MAR$table[4,]= c(GM$MAR$model12$se[c(1,2,2,3,3,4,5,6)])
GM$MAR$table[5,]= c(GM$MAR$model13$coefficients[c(1,2,2,3,4,5,6,7)])
GM$MAR$table[6,]= c(GM$MAR$model13$se[c(1,2,2,3,4,5,6,7)])
rownames(GM$MAR$table) = c("AR(2,0)","s.e.", "MAR(1,1)","s.e.", "AR(0,2)","s.e.")
colnames(GM$MAR$table) = c("int", "lag1", "lag2", "lead1", "lead2", "exo", "df", "scale" )

GM$MAR$select = selection.lag.lead(GM$cpi,x = GM$unrate,p_pseudo = 2)
selection.lag.lead(GM$cpi,x = GM$unrate,p_pseudo = 12)

GM$MAR$select5 = selection.lag.lead(GM$cpi,x = GM$unrate,p_pseudo = 5)

GM$MAR$model50 = mixed(y = GM$cpi,x = GM$unrate,p_C = 5,p_NC = 0)
GM$MAR$model41 = mixed(y = GM$cpi,x = GM$unrate,p_C = 4,p_NC = 1)
GM$MAR$model32 = mixed(y = GM$cpi,x = GM$unrate,p_C = 3,p_NC = 2)
GM$MAR$model23 = mixed(y = GM$cpi,x = GM$unrate,p_C = 2,p_NC = 3)
GM$MAR$model14 = mixed(y = GM$cpi,x = GM$unrate,p_C = 1,p_NC = 4)
GM$MAR$model05 = mixed(y = GM$cpi,x = GM$unrate,p_C = 0,p_NC = 5)

c(Box.test(GM$MAR$model50$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GM$MAR$model41$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GM$MAR$model32$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GM$MAR$model23$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GM$MAR$model14$residuals,lag = 12,type = "Ljung-Box")$p.value,
  Box.test(GM$MAR$model05$residuals,lag = 12,type = "Ljung-Box")$p.value)


GM$MAR$table2 = matrix(0,3,2)
GM$MAR$table2[,1] = GM$MAR$select$loglikelihood
GM$MAR$table2[,2] = c(Box.test(GM$MAR$model11$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GM$MAR$model12$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GM$MAR$model13$residuals,lag = 12,type = "Ljung-Box")$p.value)

GM$MAR$table3 = matrix(0,6,2)
GM$MAR$table3[,1] = GM$MAR$select5$loglikelihood
GM$MAR$table3[,2] = c(Box.test(GM$MAR$model50$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GM$MAR$model41$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GM$MAR$model32$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GM$MAR$model23$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GM$MAR$model14$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GM$MAR$model05$residuals,lag = 12,type = "Ljung-Box")$p.value)


###Great Recession Period###

n = length(GR$cpi)
nlag = 12
GR$lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        GR$cpi[1:(n - x)]))
#lags = lags[-c(1:nlag),]

GR$AIC = c()

for(i in 1:nlag){
  exo = GR$unrate[-c(1:i)]
  auxm = lm(GR$cpi[-c(1:i)]~1+GR$lags[-c(1:i),1:i]+exo)
  GR$AIC[i] = AIC(auxm)}

mixed(y = GR$cpi,x = GR$unrate,p_C = 0,p_NC = 1)
mixed(y = GR$cpi,x = GR$unrate,p_C = 1,p_NC = 0)

GR$MAR$model11 = mixed(y = GR$cpi,x = GR$unrate,p_C = 2,p_NC = 0)
GR$MAR$model12 = mixed(y = GR$cpi,x = GR$unrate,p_C = 1,p_NC = 1)
GR$MAR$model13 = mixed(y = GR$cpi,x = GR$unrate,p_C = 0,p_NC = 2)

GR$MAR$table = matrix(0,6,8)
GR$MAR$table[1,]= c(GR$MAR$model11$coefficients[c(1,2,3,4,4,5,6,7)])
GR$MAR$table[2,]= c(GR$MAR$model11$se[c(1,2,3,4,4,5,6,7)])
GR$MAR$table[3,]= c(GR$MAR$model12$coefficients[c(1,2,2,3,3,4,5,6)])
GR$MAR$table[4,]= c(GR$MAR$model12$se[c(1,2,2,3,3,4,5,6)])
GR$MAR$table[5,]= c(GR$MAR$model13$coefficients[c(1,2,2,3,4,5,6,7)])
GR$MAR$table[6,]= c(GR$MAR$model13$se[c(1,2,2,3,4,5,6,7)])
rownames(GR$MAR$table) = c("AR(2,0)","s.e.", "MAR(1,1)","s.e.", "AR(0,2)","s.e.")
colnames(GR$MAR$table) = c("int", "lag1", "lag2", "lead1", "lead2", "exo", "df", "scale" )


GR$MAR$select = selection.lag.lead(GR$cpi,x = GR$unrate,p_pseudo = 2)

GR$MAR$table2 = matrix(0,3,2)
GR$MAR$table2[,1] = GR$MAR$select$loglikelihood
GR$MAR$table2[,2] = c(Box.test(GR$MAR$model11$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GR$MAR$model12$residuals,lag = 12,type = "Ljung-Box")$p.value,
                       Box.test(GR$MAR$model13$residuals,lag = 12,type = "Ljung-Box")$p.value)


