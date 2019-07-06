#Functions


rmltest=function(data){
  #Ratio of Maximized Likelihoods (RML) test
  #the test is based on the Advances in Ranking and Selection, 
  #Multiple Comparisons, and Reliability, ed. Balakrishan, Kannan, Nagaraja. pp 65-79
  #Chapter 4. [Kundu] Discriminating between Normal and Laplace distribution
  #We are interested in testing Normality in 95 and 99 per cent interval 
  
  #Author: Igor Kindop
  #March 2017
  z95=qnorm(0.95)
  z99=qnorm(0.99)
  #minimal sample size required
  n95=47.7*z95^2
  n99=47.7*z99^2
  N=length(data)
  #Critical values for the candidate hypothesis
  #test1: H0=Normal vs. H1=Laplace
  #test2: H0=Laplace vs. H1=Normal
  
  t195=N*0.0484172-z95*sqrt(N*0.0707963);
  t199=N*0.0484172-z99*sqrt(N*0.0707963);
  t295=-N*0.0723649+z95*sqrt(N*0.25);
  t299=-N*0.0723649+z99*sqrt(N*0.25);
  Crit=as.matrix(c(t195, t199, t295,t299),1,4)
  #Test statistic for the candidate hypothesis
  Test=(log(2)-log(pi)+1)*N/2+N*(log(1/N*sum(abs(data-median(data))))-log(sd(data)));
  result=list(crit.values=Crit, size=N, nec.size=c(n95, n99), teststat=Test)
  return(result)}




#define the function we need to solve to obtain the shape parameter
f<-function(x,series) {
  (((gamma(3/x)*gamma(1/x))^0.5)/gamma(2/x)-(sd(series)/(sum(abs(series-mean(series)))/length(series))))
}


BaiTest=function(v){  
  #The function is based on the algorithm presented in the Bai(2003) 
  #and is very similar to the function used in Gaglianone and Lima.
  n=length(v)
  Wn= matrix(0,nrow=n, ncol=1)
  g_dot=matrix(0,nrow=n, ncol=3)
  
  for (k in 1:n)
  {
    g_dot[k,1] = 1
    g_dot[k,2] = -qnorm(v[k] , mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
    g_dot[k,3] = 1-(qnorm(v[k] , mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)) ^ 2
  }
    for (j in 2:n)
  {
    auxillary = 0
    for (k in 2:j)
    {
      Ck=0
      Dk=0
      for (i in k:n-1)
      {
        Ck = Ck + g_dot[i,] %*% t(g_dot[i,])*(v[i+1]-v[i])
        Dk = Dk + g_dot[i,]
      }
      
      Ck_inv = tryCatch(qr.solve(Ck,tol=1e-25),error=function(e){cat("ERROR :",c("Singularity in "), print(j), "\n")})
      auxillary = tryCatch(auxillary + t(g_dot[k,]) %*% Ck_inv %*% Dk *(v[k]-v[k-1]),error=function(e){cat("ERROR :",c("Singularity in "), print(j), "\n")})
      
      if(is.null(auxillary)==1){auxillary=NA}
    }
    
    Wn[j,1] = (n^0.5) * abs(j/n - (1/n)*auxillary)
    }
  Tn=max(Wn,na.rm=TRUE)
  result=list(teststat=Tn, testvec=Wn)
  return(result)}
