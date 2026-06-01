#Pg sampler for logistic regression with normal prior


PGNorm.logit<-function(y,Nsize,X,iter,burnin,betap,Sigmainvp){
  library(BayesLogit)
  n <-dim(X)[1]
  p <-dim(X)[2]
  #create pseudo-normal response
  kappa<-y-0.5*Nsize
  #Create of posterior mean for beta that does not change
  pbeta <-crossprod(X,kappa)
  priorFixed<-Sigmainvp%*%betap
  pbeta <- pbeta+priorFixed
  #Cholesky of prior inverse for simulating noise.
  Sigmainvchol <- chol(Sigmainvp)
  
  #Storing results.
  betavarcomp.store   <-matrix(0,iter-burnin,p+n)
  
  #Initial PG variables #assumes p = 0.5 for all y.
  myp<-sum(y)/sum(Nsize*rep(1,n)) #use for starting value
  b0   <-log(myp/(1-myp))
  omega<-rpg(n,Nsize,b0)  
  
  for(i in 1:iter){
    #generates noise X'WX (using the back solve, the correct vcov for posterior 
    #for beta is created.
    errbeta<-t(X)%*%(rnorm(n)*sqrt(omega))  + t(Sigmainvchol)%*%rnorm(p)
    bmean<- pbeta+errbeta
    Sigmainv <- crossprod(X*omega,X) + Sigmainvp
    beta<-  solve(Sigmainv,bmean) #update beta.
    Xb<-X%*%beta
    omega<-rpg(n,Nsize,Xb) 
    
    betavarcomp.store[max(1,i-burnin),]<-c(beta,omega)  
  }
  return(betavarcomp.store)  
}

