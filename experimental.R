

#next comes the experimental setup


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



