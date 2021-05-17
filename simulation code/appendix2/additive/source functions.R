################all functions
library(MASS)
library(nleqslv)
###helper
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
###generate data
data_gen <- function(n)
{
  U1 <- U2 <- U3 <- U4 <- rep(NA,n)
  for(i in 1:n)
  {
    U.temp <- mvrnorm(1,mu=cbind(0,0,0,0), Sigma=diag(c(1,1,1,1)))
    U1[i] <- U.temp[1]
    U2[i] <- U.temp[2]
    U3[i] <- U.temp[3]
    U4[i] <- U.temp[4]
  }
  U <- cbind(U1,U2,U3,U4)
  pa <- expit(U %*% aU)
  A <- rbinom(n,1,prob = pa)
  
  pm1 <- expit(m10 + m1A*A + U %*% m1U)
  pm2 <- expit(m20 + m2A*A + U %*% m2U)
  phi <- exp(tau_0+tau_A*A+U %*% tau_X)
  
  
  ##generate M1 and M2
  p11.fun <- function(phi,pm1,pm2)
  {
    if (phi==1)
    {
      return(pm1*pm2)
    }
    if(phi!=1)
    {
      s <- sqrt((1+(pm1+pm2)*(phi-1))^2+4*phi*(1-phi)*pm1*pm2)
      numer <- 1+(pm2+pm1)*(phi-1)-s
      dnom <- 2*(phi-1)
      return(numer/dnom)
    }
  }
  
  p11.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p11.vector[i] <- p11.fun(phi[i],pm1[i],pm2[i])
  }
  
  p10.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p10.vector[i] <- pm1[i]-p11.vector[i]
  }
  p01.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p01.vector[i] <- pm2[i]-p11.vector[i]
  }
  p00.vector <- rep(NA,n)
  for(i in 1:n)
  {
    p00.vector[i] <- 1-pm2[i]-pm1[i]+p11.vector[i]
  }
  
  joint <- matrix(NA, nrow = n, ncol=2)
  
  for(i in 1:n)
  {
    prob <- c(p11.vector[i],p10.vector[i],p01.vector[i],p00.vector[i])
    whichbox <- rmultinom(1,1,prob)
    index <- which(whichbox==1)
    if(index==1)
    {
      joint[i,]= c(1,1)
    }
    if(index==2)
    {
      joint[i,]= c(1,0)
    }
    if(index==3)
    {
      joint[i,]= c(0,1)
    }
    if(index==4)
    {
      joint[i,]= c(0,0)
    }
  }
  
  M1 <- joint[,1]
  M2 <- joint[,2]
  
  EY <- y0 + yA*A + yM2*M2 + yM1*M1 + yAM1*A*M1 + yAM2*A*M2 +  U %*% yU 
  
  Y <- rnorm(n,EY,sdy)
  
  Z1 <- exp(U1/2)
  Z2 <- U2/(1+exp(U1))+10
  Z3 <- (U1*U3/25+0.6)^3
  Z4 <- (U2+U4+20)^2
  
  
  dataset <- data.frame(Z1,Z2,Z3,Z4,A,M1,M2,Y,U1,U2,U3,U4)
  
  dataset
}

###joint distribution
##estimate the copular parameter
estp11<- function(fitm.mar,fitc.mar,M,C,n,newdata)
{
  fc <- expit(predict(fitc.mar,newdata = newdata))
  fm <- expit(predict(fitm.mar,newdata = newdata))
  a <- b <- c <- d <- rep(NA, n)
  
  for(i in 1:n)
  {
    if (M[i]==1 & C[i]==1)
    {
      a[i] <- 1
      b[i] <- c[i] <- d[i] <- 0
    }
    if (M[i]==0 & C[i]==1)
    {
      b[i] <- 1
      a[i] <- c[i] <- d[i] <- 0
    }
    if (M[i]==1 & C[i]==0)
    {
      c[i] <- 1
      a[i] <- b[i] <- d[i] <- 0
    }
    if (M[i]==0 & C[i]==0)
    {
      d[i] <- 1
      a[i] <- b[i] <- c[i] <- 0
    }
  }
  
  fn <- function(par)
  {
    
    f1 <- function(par)
    {
      tau0 <- par
      
      thing <- rep(NA,n)
      
      for(i in 1:n)
      {
        #phi
        phi <- exp(tau0)
        #s  
        s <- sqrt((1+(fm[i]+fc[i])*(phi-1))^2+4*(1-phi)*phi*fm[i]*fc[i])
        
        if(phi==1)
        {
          p11 <- fm[i]*fc[i]
          partialphi <- 0
        }
        if(phi!=1)
        {
          p11 <- (1+(fm[i]+fc[i])*(phi-1)-s)/(2*(phi-1))
          term11 <- fm[i]+fc[i]
          term12 <- 0.5*((1+(fm[i]+fc[i])*(phi-1))^2+4*phi*(1-phi)*fm[i]*fc[i])^(-0.5)
          term13 <- 2*(1+(fc[i]+fm[i])*(phi-1))*(fm[i]+fc[i])+4*fc[i]*fm[i]-8*fc[i]*fm[i]*phi
          term1 <- (term11-term12*term13)/(2*(phi-1))
          term2 <- p11/(phi-1)
          partialphi <- term1-term2
        }
        
        partialp11tau0 <- partialphi*phi
        
        t <- a[i]/p11-b[i]/(fc[i]-p11)-c[i]/(fm[i]-p11)+d[i]/(1-fc[i]-fm[i]+p11)
        thing[i] <- t * partialp11tau0
      }
      
      sum(thing)
    }
    
    f <- f1(par)
    f
  }
  
  p11.fun <- function(phi,pc,pm)
  {
    if (phi==1)
    {
      return(pc*pm)
    }
    if(phi!=1)
    {
      s <- sqrt((1+(pc+pm)*(phi-1))^2+4*phi*(1-phi)*pc*pm)
      numer <- 1+(pm+pc)*(phi-1)-s
      dnom <- 2*(phi-1)
      return(numer/dnom)
    }
  }
  
  p0 <- 0.5
  solution <- nleqslv(p0, fn,control = list(trace=T)) 
  
  esttau <- solution$x
  return(esttau)
}
##joint distribution of the two mediators
FMC.fun <- function(esttau,pc,pm,n)
{
  p11.fun <- function(phi,pc,pm)
  {
    if (phi==1)
    {
      return(pc*pm)
    }
    if(phi!=1)
    {
      s <- sqrt((1+(pc+pm)*(phi-1))^2+4*phi*(1-phi)*pc*pm)
      numer <- 1+(pm+pc)*(phi-1)-s
      dnom <- 2*(phi-1)
      return(numer/dnom)
    }
  }
  
  phiest <- exp(esttau)
  p11.est <- rep(NA,n)
  for(i in 1:n)
  {
    p11.est[i] <- p11.fun(phiest,pc[i],pm[i])
  }
  #print(esttau-tau_0)
  return(p11.est)
}
###original quadruply robust estimator for EIE_M1
quad_m1 <- function(dataset,fita,fity,fitm1,fitm2,copula.par)
{
  delta1.fun <- function(A,Y,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par,EYa1)
  {
    n <- length(A)
    #joint distribution
    p11.est <- FMC.fun(esttau=copula.par,pc=pm1a1,pm=pm2a1,n)
    Fm1m2 <- M1*M2*p11.est+
      (1-M1)*(1-M2)*(1-pm1a1-pm2a1+p11.est)+
      (1-M2)*M1*(pm1a1-p11.est)+
      (1-M1)*M2*(pm2a1-p11.est)
    
    fm1a1 <- dbinom(M1,1,pm1a1)
    fm1a0 <- dbinom(M1,1,pm1a0)
    fm2a1 <- dbinom(M2,1,pm2a1)
    
    da <- dbinom(A,1,proba)
    ratio <- (A/da)*((fm1a1 - fm1a0)*fm2a1/Fm1m2)
    delta1 <- ratio*Y
    
    comp1 <- ratio*(Y-EYa1)
    
    return(cbind(delta1,comp1))
  }
  
  eta1.fun <- function(fity,pm1a1,pm1a0,a,dataset)
  {
    n <- nrow(dataset)
    newdatam11 <- newdatam10 <- dataset
    newdatam11$A <- newdatam10$A <- 1
    newdatam11$M1 <- 1
    newdatam10$M1 <- 0
    pm1 <- a*pm1a1 + (1-a)*pm1a0
    eta <- predict(fity,newdata = newdatam11,type="response")*pm1 +
      predict(fity,newdata = newdatam10,type="response")*(1-pm1)
    eta
  }
  
  delta2.fun <- function(A,proba,fity,pm1a1,pm1a0,dataset)
  {
    n <- nrow(dataset)
    eta.diff <- eta1.fun(fity,pm1a1,pm1a0,a=1,dataset) - eta1.fun(fity,pm1a1,pm1a0,a=0,dataset)
    delta2 <- comp2 <- (A/proba)*eta.diff
    return(cbind(delta2,comp2))
  }
  
  eta2.fun <- function(fity,pm2a1,dataset)
  {
    n <- nrow(dataset)
    newdatam21 <- newdatam20 <- dataset
    newdatam21$A <- newdatam20$A <- 1
    newdatam21$M2 <- 1
    newdatam20$M2 <- 0
    pm2 <- pm2a1
    eta <- predict(fity,newdata = newdatam21,type="response")*pm2 +
      predict(fity,newdata = newdatam20,type="response")*(1-pm2)
    eta
  }
  
  delta3.fun <- function(A,proba,fity,pm2a1,dataset)
  {
    n <- nrow(dataset)
    da <- dbinom(A,1,proba)
    delta3 <- comp3 <- ((2*A-1)/da)*eta2.fun(fity,pm2a1,dataset)
    return(cbind(delta3,comp3))
  }
  
  gamma1.fun <- function(fity,pm1a1,pm1a0,pm2a1,a,dataset)
  {
    n <- nrow(dataset)
    newdata11 <- newdata10 <- newdata01 <- newdata00 <- dataset
    newdata11$A <- newdata10$A <- newdata01$A <- newdata00$A <- 1
    newdata11$M1 <- newdata10$M1 <- newdata11$M2 <- newdata01$M2 <- 1
    newdata10$M2 <- newdata00$M2 <- newdata01$M1 <- newdata00$M1 <- 0 
    
    pm1 <- a*pm1a1 + (1-a)*pm1a0
    pm2 <- pm2a1
    
    gamma <- predict(fity,newdata = newdata11,type="response")*pm1*pm2 +
      predict(fity,newdata = newdata10,type="response")*pm1*(1-pm2) +
      predict(fity,newdata = newdata01,type="response")*(1-pm1)*pm2 +
      predict(fity,newdata = newdata00,type="response")*(1-pm1)*(1-pm2) 
    
    gamma
  }
  
  delta4.fun <- function(A,proba,fity,pm1a1,pm1a0,pm2a1,dataset)
  {
    n <- nrow(dataset)
    da <- dbinom(A,1,proba)
    gamma1 <- gamma1.fun(fity,pm1a1,pm1a0,pm2a1,1,dataset)
    gamma0 <- gamma1.fun(fity,pm1a1,pm1a0,pm2a1,0,dataset)
    delta4 <-  gamma1 - gamma0
    comp4 <- delta4 - (2*A/da)*gamma1 + (1/da)*gamma0
    
    return(cbind(delta4,comp4))
  }
  
  n <- nrow(dataset)
  #initiate parameters
  A <- dataset$A
  M1 <- dataset$M1
  M2 <- dataset$M2
  Y <- dataset$Y
  #  dataset$R1 <- 0
  
  proba <- fitted(fita)
  newdata1 <- newdata0 <- dataset
  newdata1$A <- 1
  newdata0$A <- 0
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  pm1a0 <- expit(predict(fitm1,newdata = newdata0))
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  EYa1 <- predict(fity,newdata = newdata1,type="response")
  
  one <- delta1.fun(A,Y,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par,EYa1)
  two <- delta2.fun(A,proba,fity,pm1a1,pm1a0,newdata1)
  three <- delta3.fun(A,proba,fity,pm2a1,newdata1)
  four <- delta4.fun(A,proba,fity,pm1a1,pm1a0,pm2a1,newdata1)
  
  delta1 <- mean(one[,1])
  delta2 <- mean(two[,1])
  delta3 <- mean(three[,1])
  delta4 <- mean(four[,1])
  
  delta.quad <- mean(one[,2] + two[,2] + three[,2] + four[,2])
  
  c(delta1,delta2,delta3,delta4,delta.quad)
}
###stablized quadruply robust estimator for EIE_M1
quad_m1_stable <- function(dataset,fita,fity,fitm1,fitm2,copula.par)
{
  delta1.fun <- function(A,Y,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par,EYa1)
  {
    n <- length(A)
    #joint distribution
    p11.est <- FMC.fun(esttau=copula.par,pc=pm1a1,pm=pm2a1,n)
    Fm1m2 <- M1*M2*p11.est+
      (1-M1)*(1-M2)*(1-pm1a1-pm2a1+p11.est)+
      (1-M2)*M1*(pm1a1-p11.est)+
      (1-M1)*M2*(pm2a1-p11.est)
    
    fm1a1 <- dbinom(M1,1,pm1a1)
    fm1a0 <- dbinom(M1,1,pm1a0)
    fm2a1 <- dbinom(M2,1,pm2a1)
    
    da <- dbinom(A,1,proba)
    ratio <- (A/da)*((fm1a1 - fm1a0)*fm2a1/Fm1m2)
    delta1 <- ratio*Y
    
    comp1 <- ratio*(Y-EYa1)
    
    return(cbind(delta1,comp1))
  }
  
  eta1.fun <- function(fity,proba,pm1a1,pm1a0,pm2a1,copula.par,newdata1,a)
  {
    n <- nrow(newdata1)
    newdatam11 <- newdatam10 <- newdata1
    #  newdatam11$A <- newdatam10$A <- 1
    newdatam11$M1 <- 1
    newdatam10$M1 <- 0
    
    newdatam11$R1 <- ratio.fun1_R1(newdatam11$M1,newdatam11$M2,pm1a1,pm1a0,pm2a1,copula.par)
    newdatam10$R1 <- ratio.fun1_R1(newdatam10$M1,newdatam10$M2,pm1a1,pm1a0,pm2a1,copula.par)
    
    pm1 <- a*pm1a1 + (1-a)*pm1a0
    eta <- predict(fity,newdata = newdatam11,type="response")*pm1 +
      predict(fity,newdata = newdatam10,type="response")*(1-pm1)
    eta
  }
  
  delta2.fun <- function(A,proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
  {
    n <- nrow(newdata1)
    eta.diff <- eta1.fun(fity,proba,pm1a1,pm1a0,pm2a1,copula.par,newdata1,a=1) - eta1.fun(fity,proba,pm1a1,pm1a0,pm2a1,copula.par,newdata1,a=0)
    delta2 <- comp2 <- (A/proba)*eta.diff
    return(cbind(delta2,comp2))
  }
  
  eta2.fun <- function(proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
  {
    n <- nrow(newdata1)
    newdatam21 <- newdatam20 <- newdata1
    #  newdatam21$A <- newdatam20$A <- 1
    newdatam21$M2 <- 1
    newdatam20$M2 <- 0
    
    newdatam21$R1 <- ratio.fun1_R1(newdatam21$M1,newdatam21$M2,pm1a1,pm1a0,pm2a1,copula.par)
    newdatam20$R1 <- ratio.fun1_R1(newdatam20$M1,newdatam20$M2,pm1a1,pm1a0,pm2a1,copula.par)
    
    pm2 <- pm2a1
    eta <- predict(fity,newdata = newdatam21,type="response")*pm2 +
      predict(fity,newdata = newdatam20,type="response")*(1-pm2)
    eta
  }
  
  delta3.fun <- function(A,proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
  {
    n <- nrow(newdata1)
    da <- dbinom(A,1,proba)
    delta3 <- comp3 <- ((2*A-1)/da)*eta2.fun(proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
    return(cbind(delta3,comp3))
  }
  
  gamma1.fun <- function(proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par,a)
  {
    n <- nrow(newdata1)
    newdata11 <- newdata10 <- newdata01 <- newdata00 <- newdata1
    #  newdata11$A <- newdata10$A <- newdata01$A <- newdata00$A <- 1
    newdata11$M1 <- newdata10$M1 <- newdata11$M2 <- newdata01$M2 <- 1
    newdata10$M2 <- newdata00$M2 <- newdata01$M1 <- newdata00$M1 <- 0 
    
    pm1 <- a*pm1a1 + (1-a)*pm1a0
    pm2 <- pm2a1
    
    newdata11$R1 <- ratio.fun1_R1(newdata11$M1,newdata11$M2,pm1a1,pm1a0,pm2a1,copula.par)
    newdata10$R1 <- ratio.fun1_R1(newdata10$M1,newdata10$M2,pm1a1,pm1a0,pm2a1,copula.par)
    newdata01$R1 <- ratio.fun1_R1(newdata01$M1,newdata01$M2,pm1a1,pm1a0,pm2a1,copula.par)
    newdata00$R1 <- ratio.fun1_R1(newdata00$M1,newdata00$M2,pm1a1,pm1a0,pm2a1,copula.par)
    
    gamma <- predict(fity,newdata = newdata11,type="response")*pm1*pm2 +
      predict(fity,newdata = newdata10,type="response")*pm1*(1-pm2) +
      predict(fity,newdata = newdata01,type="response")*(1-pm1)*pm2 +
      predict(fity,newdata = newdata00,type="response")*(1-pm1)*(1-pm2) 
    
    gamma
  }
  
  delta4.fun <- function(A,proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
  {
    n <- nrow(newdata1)
    da <- dbinom(A,1,proba)
    gamma1 <- gamma1.fun(proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par,a=1)
    gamma0 <- gamma1.fun(proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par,a=0)
    delta4 <-  gamma1 - gamma0
    comp4 <- delta4 - (2*A/da)*gamma1 + (1/da)*gamma0
    
    return(cbind(delta4,comp4))
  }
  
  n <- nrow(dataset)
  #initiate parameters
  A <- dataset$A
  M1 <- dataset$M1
  M2 <- dataset$M2
  Y <- dataset$Y
  #  dataset$R1 <- 0
  
  proba <- fitted(fita)
  newdata1 <- newdata0 <- dataset
  newdata1$A <- 1
  newdata0$A <- 0
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  pm1a0 <- expit(predict(fitm1,newdata = newdata0))
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  
  newdata1$R1 <- ratio.fun1_R1(M1,M2,pm1a1,pm1a0,pm2a1,copula.par)
  #  newdata0$R0 <- ratio.fun1(0,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par)
  
  EYa1 <- predict(fity,newdata = newdata1,type="response")
  
  one <- delta1.fun(A,Y,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par,EYa1)
  two <- delta2.fun(A,proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
  three <- delta3.fun(A,proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
  four <- delta4.fun(A,proba,fity,pm1a1,pm1a0,pm2a1,newdata1,copula.par)
  
  delta1 <- mean(one[,1],na.rm = T)
  delta2 <- mean(two[,1],na.rm = T)
  delta3 <- mean(three[,1],na.rm = T)
  delta4 <- mean(four[,1],na.rm = T)
  
  delta.quad <- mean(one[,2] + two[,2] + three[,2] + four[,2],na.rm = T)
  
  c(delta1,delta2,delta3,delta4,delta.quad)
}

ratio.fun1_R1 <- function(M1,M2,pm1a1,pm1a0,pm2a1,copula.par)
{
  n <- length(M1)
  #joint distribution
  p11.est <- FMC.fun(esttau=copula.par,pc=pm1a1,pm=pm2a1,n)
  Fm1m2 <- M1*M2*p11.est+
    (1-M1)*(1-M2)*(1-pm1a1-pm2a1+p11.est)+
    (1-M2)*M1*(pm1a1-p11.est)+
    (1-M1)*M2*(pm2a1-p11.est)
  
  fm1a1 <- dbinom(M1,1,pm1a1)
  fm1a0 <- dbinom(M1,1,pm1a0)
  fm2a1 <- dbinom(M2,1,pm2a1)
  
  ratio <- ((fm1a1 - fm1a0)*fm2a1/Fm1m2)
  
  return(ratio)
}

stable1.m1 <- function(dataset,fita,M1T=T,M2T=T,YT=T)
{
  n <- nrow(dataset)
  A <- dataset$A
  M1 <- dataset$M1
  M2 <- dataset$M2
  varname.m1 <- ifelse(M1T,"U","Z")
  varname.m2 <- ifelse(M2T,"U","Z")
  varname.y <- ifelse(YT,"U","Z")
  for.m1 <- paste0("cbind(M1,1-M1)~",paste(c("A",paste0(rep(varname.m1,4),1:4)),collapse = "+"))
  for.m2 <- paste0("cbind(M2,1-M2)~",paste(c("A",paste0(rep(varname.m2,4),1:4)),collapse = "+"))
  for.y <- paste0("Y~",paste(c("A","M1","M2","R1","A*M1","A*M2",paste0(rep(varname.y,4),1:4)),collapse = "+"))
  
  proba <- fitted(fita)
  dataset$weighta <- (A*proba+(1-A)*(1-proba))^(-1)
  fitm1 <- glm(for.m1,family = binomial(link = "logit"),data = dataset,
               weights = weighta)
  fitm2 <- glm(for.m2,family = binomial(link = "logit"),data = dataset,
               weights = weighta)
  
  esttau <- estp11(fitm2,fitm1,M2,M1,nrow(dataset),dataset)
  
  newdata1 <- newdata0 <- dataset
  newdata1$A <- 1
  newdata0$A <- 0
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  pm1a0 <- expit(predict(fitm1,newdata = newdata0))
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  
  R1 <- ratio.fun1_R1(M1,M2,pm1a1,pm1a0,pm2a1,copula.par=esttau)
  dataset$R1 <- as.numeric(R1)
  
  fity <- lm(for.y,data = dataset,weights=weighta)
  quad.stable <- quad_m1_stable(dataset,fita,fity,fitm1,fitm2,esttau)
  
  return(quad.stable[5])
}
