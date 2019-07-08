a = selection.lag(y = GGI$cpi,x = GGI$unrate,p_max = 12)$aic
b = VARselect(y = GGI$cpi,lag.max = 12,exogen = GGI$unrate,type = )$criteria[1,]

plot(1:13,c(a), type = "l")
lines(c(b[1],b)+3,col = "red")


aicA = c()
aicB = c()
for (i in 1:15){
  a = VARselect(y = cbind(GGI$cpi, GGI$unrate),lag.max = i,type = "const")$criteria[1,]
  aicA[i] = which.min(a)
  aicB[i] = a[aicA[i]]
  }

###tests on Normality of the residuals

jarque.bera.test(GGI$MAR$model12$residuals)

rmltest(GGI$causalmodel$residuals) #reject Normality and do not reject Gaussian
rmltest(GGM$causalmodel$residuals)
rmltest(GGR$causalmodel$residuals)


par = 5

shape1 = optim(par = par,fn = (f), series = GGI$causalmodel$residuals)
shape2 =optim(par = par,fn = (f), series = GGM$causalmodel$residuals)
shape3 =optim(par = par,fn = (f), series = GGR$causalmodel$residuals)


stdresA = (GGI$causalmodel$residuals - mean(GGI$causalmodel$residuals))/sd(GGI$causalmodel$residuals)
stdresA = sort(pnorm(stdresA))
BaiTest(stdresA) #1.921316

#1.94/ 2.22/ 2.80 for 10%/5%/1% significance levels

stdresB = (GGM$causalmodel$residuals - mean(GGM$causalmodel$residuals))/sd(GGM$causalmodel$residuals)
stdresB = sort(pnorm(stdresB))
BaiTest(stdresB) #3.234138

stdresC = (GGR$causalmodel$residuals - mean(GGR$causalmodel$residuals))/sd(GGR$causalmodel$residuals)
stdresC = sort(pnorm(stdresC))
BaiTest(stdresC) #2.013473

####Here I estimate the shape parameter of the EPD

shapeA = uniroot(f2, c(0.1,6), series=GGI$causalmodel$residuals)
shapeB = uniroot(f2, c(0.1,6), series=GGM$causalmodel$residuals)
shapeC = uniroot(f2, c(0.1,6), series=GGR$causalmodel$residuals)

#wild bootstrap to validate our estimates of the shape parameter

btA = c()
btB = c()
btC = c()
btAo = c()
btBo = c()
btCo = c()
bt0 = c()

n1 = length(GGI$causalmodel$residuals)
n2 = length(GGM$causalmodel$residuals)
n3 = length(GGR$causalmodel$residuals)
par = 5

for (i in 1:1000){
  x = rnorm(length(GGI$cpi)) #variance of the modifier is one
  bt0[i] = uniroot(f2, c(0.1,6), series=GGI$cpi*x)$root
  x = rnorm(n1) #variance of the modifier is one
  btA[i] = uniroot(f2, c(0.1,6), series=GGI$causalmodel$residuals*x)$root
  btAo[i] = optim(par = par,fn = (f), series = GGI$causalmodel$residuals*x)$par
  
  x = rnorm(n2)
  btB[i] = uniroot(f2, c(0.1,6), series=GGM$causalmodel$residuals*x)$root
  btBo[i] = optim(par = par,fn = (f), series = GGM$causalmodel$residuals*x)$par
  
  x = rnorm(n3)
  btC[i] = uniroot(f2, c(0.1,6), series=GGR$causalmodel$residuals*x)$root
  btCo[i] = optim(par = par,fn = (f), series = GGR$causalmodel$residuals*x)$par
  
  }

###simple nonparametric bootstrap and the EPDs shape

library(boot)
z = function(d,w) uniroot(f2, c(0.1,6), d[w])$root
resultsA = boot(data = GGI$causalmodel$residuals, 
               statistic = z, R = 1000, stype = "i")
boot.ci(boot.out = resultsA, conf = 0.95, type=c("norm"))

resultsB = boot(data = GGM$causalmodel$residuals, 
                statistic = z, R = 1000, stype = "i")
boot.ci(boot.out = resultsB, conf = 0.95, type=c("norm"))

resultsC = boot(data = GGR$causalmodel$residuals, 
                statistic = z, R = 1000, stype = "i")
boot.ci(boot.out = resultsC, conf = 0.95, type=c("norm"))


