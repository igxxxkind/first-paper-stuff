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

