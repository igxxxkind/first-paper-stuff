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
# LAD=function(p, yt, r, order){
#   #order - causal, noncausal roots
#   n = length(yt)
#   nlag = order[1]
#   nlead = order[2]
#   lags = c()
#   leads = c()
#   if(nlag>0){
#     lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
#                                         yt[1:(n - x)]))}else{lags = 1}
#   if(nlead>0){
#     leads = sapply(1:nlead, function(x) c(yt[-c(1:x)],
#                                           rep(NA, length.out = x)))}else{leads = 1}
#   if(nlag==0&&nlead==0){
#     regs = rep(1,length.out = n)}else if(nlead==0){
#       regs = cbind(leads, lags)}else if (nlag==0){
#         regs = cbind(lags, leads)}else {
#           regs = cbind(1,lags,leads)}
#   
#   if(!is.null(r)){
#     regs = cbind(regs,r)}
#   if(length(p)!=dim(regs)[2]){p = p[-length(p)]}
#   rhs = p%*%t(regs)
#   
#   cond=yt[!is.na(rhs)]-rhs[!is.na(rhs)]
#   stat=sum(abs(cond))
#   return(stat)
# }



LAD=function(par, yt, r,order){
  #order - (causal, noncausal) roots
  n = length(yt)
  nlag = order[1]
  nlead = order[2]
  lags = c()
  leads = c()
  r = matrix(r,n,)
  
  if(nlag>0){
    lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        yt[1:(n - x)]))}else{lags = 1}
  if(nlead>0){
    leads = sapply(1:nlead, function(x) c(yt[-c(1:x)],
                                          rep(NA, length.out = x)))}else{leads = 1}
  if(nlag==0&&nlead==0){
    regs = rep(1,length.out = n)}else if(nlead==0){
      regs = cbind(leads, lags)}else if (nlag==0){
        regs = cbind(lags, leads)}else {
          #i have to think about this cases later
          regs = cbind(lags[,dim(lags)[2]:1],yt,leads)}
  #regs[lag2,lag1,lag0,lead1,lead2]
  
  nc = c(1,-par[1+c(1:nlead)]) #noncausal lead component
  cc = c(1,-par[length(nc)+(1:nlag)]) #causal lag component
  
  marregs = matrix(0,n,1+nlead)
      for(i in 1:(nlag+1)){
      marregs = marregs+cc[i]*regs[,dim(regs)[2]-(nlead:0)-i+1]
  }
  lhs = nc%*%t(marregs)  
  
  if(!is.null(r)){
  cond=lhs[!is.na(lhs)]-par[1] - 
    par[c(1+nlag+nlead+1):length(par)]%*%t(r[!is.na(lhs),])}else{
  cond=lhs[!is.na(lhs)]-par[1]}  
  stat=sum(abs(cond))
  return(stat)
}

LADest = function(yt, exoreg, order){
  exoreg = matrix(exoreg,length(yt),)
  p = rep(0,1+sum(order)+dim(exoreg)[2])
  
  check5 = optim(par = p,fn = LAD,
                 yt = yt, 
                 r = exoreg,
                 order = order,
                 control = list(maxit = 2000),hessian = TRUE)
  print(paste("Convergence report:",check5$convergence,sep = " "))
  coefficients = as.data.frame(t(check5$par))
  LogLike = check5$value
  se = sqrt(diag(solve(check5$hessian)))
  colnames(coefficients)[1]= c("const")
  colnames(coefficients)[1+1:order[2]] = paste("Lead",1:order[2], sep = "")
  colnames(coefficients)[1+order[2]+1:order[1]] = paste("Lag",1:order[1], sep = "")
  colnames(coefficients)[-c(1:(sum(order)+1))] = paste("Exo",1:dim(exoreg)[2],sep = "")
  return(list(coefficients = coefficients,
              LogLike = LogLike,
              s.e. = se))
  }


TMLE=function(p, yt, r, order){
  pn = length(p)
  sigma=p[pn-1]
  lambda=p[pn]
  par = p[-c(pn-1,pn)]
  #order - (causal, noncausal) roots
  n = length(yt)
  nlag = order[1]
  nlead = order[2]
  lags = c()
  leads = c()
  r = matrix(r,n,)
  
  if(nlag>0){
    lags = sapply(1:nlag, function(x) c(rep(NA, length.out = x), 
                                        yt[1:(n - x)]))}else{lags = 1}
  if(nlead>0){
    leads = sapply(1:nlead, function(x) c(yt[-c(1:x)],
                                          rep(NA, length.out = x)))}else{leads = 1}
  if(nlag==0&&nlead==0){
    regs = rep(1,length.out = n)}else if(nlead==0){
      regs = cbind(leads, lags)}else if (nlag==0){
        regs = cbind(lags, leads)}else {
          #i have to think about this cases later
          regs = cbind(lags[,dim(lags)[2]:1],yt,leads)}
  #regs[lag2,lag1,lag0,lead1,lead2]
  
  nc = c(1,-par[1+c(1:nlead)]) #noncausal lead component
  cc = c(1,-par[length(nc)+(1:nlag)]) #causal lag component
  
  marregs = matrix(0,n,1+nlead)
  for(i in 1:(nlag+1)){
    marregs = marregs+cc[i]*regs[,dim(regs)[2]-(nlead:0)-i+1]
  }
  lhs = nc%*%t(marregs)  
  
  if(!is.null(r)){
    cond=lhs[!is.na(lhs)]-par[1] - 
      par[c(1+nlag+nlead+1):length(par)]%*%t(r[!is.na(lhs),])}else{
        cond=lhs[!is.na(lhs)]-par[1]}  
  
  ll1=log(gamma((lambda+1)/2)/gamma(lambda/2)/sqrt(pi))-0.5*log(lambda-2)-log(sigma)
  ll2=-(lambda+1)/2*log(1+cond^2/sigma^2/(lambda-2))
  loglike=-(sum(ll1+ll2))
  return(loglike)
}

TMLEest = function(yt, exoreg, order){
  exoreg = matrix(exoreg,length(yt),)
  p = rep(0,1+sum(order)+dim(exoreg)[2])
  p = c(p,3,3)
  
  check5 = optim(par = p,fn = TMLE,
                 yt = yt, 
                 r = exoreg,
                 order = order,
                 control = list(maxit = 5000),hessian = TRUE)
  print(paste("Convergence report:",check5$convergence,sep = " "))
  n = length(check5$par)
  coefficients = as.data.frame(t(check5$par))
  LogLike = check5$value
  se = sqrt(diag(solve(check5$hessian)))
  colnames(coefficients)[1]= c("const")
  colnames(coefficients)[1+1:order[2]] = paste("Lead",1:order[2], sep = "")
  colnames(coefficients)[1+order[2]+1:order[1]] = paste("Lag",1:order[1], sep = "")
  colnames(coefficients)[-c(1:(sum(order)+1),n-1,n)] = paste("Exo",1:dim(exoreg)[2],sep = "")
  colnames(coefficients)[c(n-1,n)] = c("sigma", "lambda")
  return(list(coefficients = coefficients,
              LogLike = LogLike,
              s.e. = se))
}

