# cpi = GI$cpi
# unrate = GI$unrate
# range = 20
# order = 30

diagnostics = function(cpi, unrate, range, order) {
  n = length(cpi)
  lags = list()
  bgtestbattery = list()
  lags[[1]] = cbind(cpi[-(1)],cpi[-(n)])
  lmaux = list()
  lmaux[[1]] = lm(lags[[1]][,1]~1+lags[[1]][,-1]+unrate[-1])
  bgtestbattery[[1]] = bgtest(formula = lmaux[[1]],order = order)
  bgtest = bgtestbattery[[1]]$statistic
  bgtestpval = bgtestbattery[[1]]$p.value
  lbtest = Box.test(x = lmaux[[1]]$residuals,lag = order,type = "Ljung-Box")$statistic
  lbtestpval = Box.test(x = lmaux[[1]]$residuals,lag = order,type = "Ljung-Box")$p.value
  jbtest = jarque.bera.test(x = lmaux[[1]]$residuals)$statistic
  jbtestpval = jarque.bera.test(x = lmaux[[1]]$residuals)$p.value
  
  for(i in 1:range){
    #there may be a problem with the bgtest
    lags[[i+1]] = cbind(cpi[-(1:(i+1))],lags[[i]][-dim(lags[[i]])[1],])
    lmaux[[i+1]] = lm(lags[[i+1]][,1]~1+lags[[i+1]][,-1]+unrate[-(1:(i+1))])
    bgtestbattery[[i+1]] = bgtest(formula = lmaux[[i+1]],order = order,fill = NA)
    bgtest[i+1] = bgtestbattery[[i+1]]$statistic
    bgtestpval[i+1] = bgtestbattery[[i+1]]$p.value
    lbtest[i+1] = Box.test(x = lmaux[[i+1]]$residuals,lag = order,type = "Ljung-Box")$statistic
    lbtestpval[i+1] = Box.test(x = lmaux[[i+1]]$residuals,lag = order,type = "Ljung-Box")$p.value
    jbtest[i+1] = jarque.bera.test(x = lmaux[[i+1]]$residuals)$statistic
    jbtestpval[i+1] = jarque.bera.test(x = lmaux[[i+1]]$residuals)$p.value
    }
  result = list(JBstat = jbtest, JBp.value = jbtestpval, 
                LBtest = lbtest, LBp.value = lbtestpval,
                BGtest = bgtest,BGp.value = bgtestpval,
                lags = lags, LSoutput = lmaux, BGtestoutput = bgtestbattery)
  return(result)
}




modelfit2 = function(data,l,exo = 1){
  #the functions picks the data and chooses the most appropriate lag structure
  #according to AIC or BIC wthin the given number of lags
  
  #Written by: Igor Kindop
  #20/04/2017 
  
  m=length(data)
  lags=matrix(1,m,1)
  if(length(exo)!=1){lags = cbind(lags,exo)}
  aic=matrix(0,l,1)
  bic=aic
  phi = list()
  for(i in 1:l)
  { lags=cbind(lags[-1,],data[-((m-i+1):m)])
  phi[[i]]=lm(data[(i+1):m]~lags[,1:min(dim(lags))])
  aic[i,]=AIC(phi[[i]])
  bic[i,]=BIC(phi[[i]])
  }
  optaic=which.min(aic)
  optbic=which.min(bic)
  return(list(aic = aic, bic = bic, optAIC=optaic, optBIC=optbic))
}


# 

