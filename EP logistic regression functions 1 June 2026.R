#EP example to improve logit-normal paper.
#1 June 2026


#Function 1: Numerical Integration.
EPNI.logit<-function(response,n,X,iter,epsilon,Prior.mean,Prior.Sinv){
  Ny<-length(response) #size of dataset.
  p<-dim(X)[2]
  
  
  #Natural parameters for Normal prior without hyper-parameters.
  PriorSinvmu <- Prior.Sinv%*%Prior.mean
  
  
  Sigmainvmu <-matrix(0,p,Ny+1) #natural parameter \Sigma^{-1}mu
  Sigmainvmu[,Ny+1]<-PriorSinvmu
  Sigmainv  <-rep(list(diag(p)),Ny)  #natural parameter \Sigma^{-1}, stored as list
  Sigmainv[[Ny+1]]<-Prior.Sinv
  
  #Previous parameters of g.
  Sigmainvmu0<-rowSums(Sigmainvmu)
  Sigmainv0   <-Reduce("+",Sigmainv)
  
  #g_0(\bm \beta) need not be updated as normal with no hyper-parameters. 
  #Note this includes the flat prior case.
  
  #function for tilted distribution
  tilt.dist      <- function(x){
    gnoti          <-dnorm(x,mean=Mnoti,sd=sqrt(Vnoti))
    pr              <-(1+exp(-x))^(-1)
    like           <-dbinom(response[i],n[i],pr)
    result <-gnoti*like
    return(result)
  }
  
  #loop for updating g_i(\bm \beta) i= 1, ..., n.
  for(j in 1:iter){
    for(i in 1:Ny){
      Sigmainvmunoti <-rowSums(Sigmainvmu)  - Sigmainvmu[,i] #Natural parameter \Sigma_{-i}^{-1}\mu_{-i}
      Sigmainvnoti   <-Reduce("+",Sigmainv) - Sigmainv[[i]]  #Natural parameter \Sigma_{-i}^{-1}\mu_{-i}
      Sigmanoti      <-solve(Sigmainvnoti)                   #parameter \Sigma_{-i}
      munoti         <-Sigmanoti%*%Sigmainvmunoti            #parameter \mu_{-i}
      Mnoti          <-t(X[i,])%*%munoti                     #M_{-i}
      Vnoti          <-t(X[i,])%*%Sigmanoti%*%X[i,]          #V_{-i}
      
      Mnoti  <-as.numeric(Mnoti)
      Vnoti  <-as.numeric(Vnoti)
      
      #Moment matching.
      E0<-integrate(f= function(x) { tilt.dist(x)}, lower=Mnoti-10*sqrt(Vnoti), upper=Mnoti+10*sqrt(Vnoti))
      E1<-integrate(f= function(x) {x*tilt.dist(x)}, lower=Mnoti-10*sqrt(Vnoti), upper=Mnoti+10*sqrt(Vnoti))
      E2<-integrate(f= function(x) {x*x*tilt.dist(x)}, lower=Mnoti-10*sqrt(Vnoti), upper=Mnoti+10*sqrt(Vnoti))  
      M <- E1$value/E0$value      
      V <- E2$value/E0$value - M*M
      
      #Update g_i
      MiViinv <- M/V - Mnoti/Vnoti
      Viinv   <- 1/V - 1/Vnoti
      
      #transform back to beta scale.
      Sigmainvmu[,i] <-X[i,]*MiViinv  #natural parameter \Sigma^{-1}mu
      Sigmainv[[i]]  <-X[i,]%*%t(X[i,])*Viinv
    }
    #Note by the way the previous lines of code have been written, step six has been implicitly. 
    
    #Checking whether to stop iterations.
    currentSinvmu <-rowSums(Sigmainvmu)   
    currentSinv   <-Reduce("+",Sigmainv)
    
    diff1  <- sqrt((currentSinvmu-Sigmainvmu0)^2)/(abs(currentSinvmu)+0.01)
    diff2  <- sqrt((currentSinv-Sigmainv0)^2)/(abs(currentSinv)+0.01)
    diff.all<-c(diff1,diff2)
    if(max(diff.all) < epsilon) break else Sigmainvmu0 <- currentSinvmu; Sigmainv0 <- currentSinv 
  }
  
  #Final mean and variance-covariance matrix of g(\beta)
  Sigma <-solve(currentSinv)
  mu    <-Sigma%*%currentSinvmu
  
  
  #Storing and returning results.
  param<-list(mu,Sigma,j)
  names(param)<-c('betahat','Sigma','iter_break')
  return(param)  
}

#Function 1: Numerical Integration.
#code written assuming data is binary, not binomial.

#This functions depend on Rcpp function
library(Rcpp)
sourceCpp('logitnormalmomentsEPY1function.cpp')

EPLN.binlogit<-function(response,X,iter,epsilon,Prior.mean,Prior.Sinv){
  Ny<-length(response) #size of dataset.
  p<-dim(X)[2]
  #sign indicator needing for moment matching.
  mysigns<-2*response -1 
  
  #Natural parameters for Normal prior without hyper-parameters.
  PriorSinvmu <- Prior.Sinv%*%Prior.mean
  
  Sigmainvmu <-matrix(0,p,Ny+1) #natural parameter \Sigma^{-1}mu
  Sigmainvmu[,Ny+1]<-PriorSinvmu
  Sigmainv  <-rep(list(diag(p)),Ny)  #natural parameter \Sigma^{-1}, stored as list
  Sigmainv[[Ny+1]]<-Prior.Sinv
  
  #Previous parameters of g.
  Sigmainvmu0<-rowSums(Sigmainvmu)
  Sigmainv0   <-Reduce("+",Sigmainv)
  
  #g_0(\bm \beta) need not be updated as normal with no hyper-parameters. Note this included
  #the flat prior case.
  
  
  #loop for updating g_i(\bm \beta) i= 1, ..., Ny.
  for(j in 1:iter){
    for(i in 1:Ny){
      Sigmainvmunoti <-rowSums(Sigmainvmu)  - Sigmainvmu[,i] #Natural parameter \Sigma_{-i}^{-1}\mu_{-i}
      Sigmainvnoti   <-Reduce("+",Sigmainv) - Sigmainv[[i]]  #Natural parameter \Sigma_{-i}^{-1}\mu_{-i}
      Sigmanoti      <-solve(Sigmainvnoti)                   #parameter \Sigma_{-i}
      munoti         <-Sigmanoti%*%Sigmainvmunoti            #parameter \mu_{-i}
      Mnoti          <-t(X[i,])%*%munoti                     #M_{-i}
      Vnoti          <-t(X[i,])%*%Sigmanoti%*%X[i,]          #V_{-i}
      
      Mnoti  <-as.numeric(Mnoti)
      Vnoti  <-as.numeric(Vnoti)
      
      #Moment matching.
      E012<-EPint_logitnormal(mu=Mnoti,sigma=sqrt(Vnoti),signY=mysigns[i])
      M <- E012[2]/E012[1]     
      V <- E012[3]/E012[1] - M*M
      
      #Update g_i
      MiViinv <- M/V - Mnoti/Vnoti
      Viinv   <- 1/V - 1/Vnoti
      
      #transform back to beta scale.
      Sigmainvmu[,i] <-X[i,]*MiViinv  #natural parameter \Sigma^{-1}mu
      Sigmainv[[i]]  <-X[i,]%*%t(X[i,])*Viinv
    }
    #Note by the way the previous lines of code have been written, step six has been implicitly. 
    
    #Checking whether to stop iterations.
    currentSinvmu <-rowSums(Sigmainvmu)   
    currentSinv   <-Reduce("+",Sigmainv)
    
    diff1  <- sqrt((currentSinvmu-Sigmainvmu0)^2)/(abs(currentSinvmu)+0.01)
    diff2  <- sqrt((currentSinv-Sigmainv0)^2)/(abs(currentSinv)+0.01)
    diff.all<-c(diff1,diff2)
    if(max(diff.all) < epsilon) break else Sigmainvmu0 <- currentSinvmu; Sigmainv0 <- currentSinv 
  }
  
  #Final mean and variance-covariance matrix of g(\beta)
  Sigma <-solve(currentSinv)
  mu    <-Sigma%*%currentSinvmu
  
  
  #Storing and returning results.
  param<-list(mu,Sigma,j)
  names(param)<-c('betahat','Sigma','iter_break')
  return(param)  
}