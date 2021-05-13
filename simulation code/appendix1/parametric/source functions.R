
################all functions
library(VineCopula)
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
  
  EY <- y0 + yA*A + yM2*M2 + yM1*M1 + U %*% yU 
  
  Y <- rnorm(n,EY,sdy)
  
  Z1 <- exp(U1/2)
  Z2 <- U2/(1+exp(U1))+10
  Z3 <- (U1*U3/25+0.6)^3
  Z4 <- (U2+U4+20)^2
  
  
  dataset <- data.frame(Z1,Z2,Z3,Z4,A,M1,M2,Y,U1,U2,U3,U4)
  
  dataset
}

###joint distribution
estp11<- function(fitm.mar,fitc.mar,M,C,n,newdata)
{
  #psuedo MLE
  
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

###NIE
trip_eta.fun <- function(pm1a1,pm1a0,pm2a1,pm2a0,fity,dataset,e,estar,copula.par)
{
  n <- nrow(dataset)
  pm1 <- estar*pm1a1 + (1 - estar)*pm1a0
  pm2 <- estar*pm2a1 + (1 - estar)*pm2a0
  
  p11.est <- FMC.fun(esttau=copula.par,pc=pm1,pm=pm2,n)
  data11 <- data10 <- data01 <- data00 <- dataset
  data11$A <- data10$A <- data01$A <- data00$A <- e
  data11$M1 <- data11$M2 <- data10$M1 <- data01$M2 <- 1
  data00$M1 <- data00$M2 <- data10$M2 <- data01$M1 <- 0
  
  
  eta <-   predict(fity,newdata = data11)*p11.est +
    predict(fity,newdata = data00)*(1-pm1-pm2+p11.est) +
    predict(fity,newdata = data10)*(pm1-p11.est) +
    predict(fity,newdata = data01)*(pm2-p11.est)
  
  eta
}

med.distrn <- function(pm1a1,pm1a0,pm2a1,pm2a0,estar,M1,M2,copula.par)
{
  n <- length(M1)
  pm1 <- estar*pm1a1 + (1 - estar)*pm1a0
  pm2 <- estar*pm2a1 + (1 - estar)*pm2a0
  
  p11.est <- FMC.fun(esttau=copula.par,pc=pm1,pm=pm2,n)
  Fm1m2 <- M1*M2*p11.est+
    (1-M1)*(1-M2)*(1-pm1-pm2+p11.est)+
    (1-M2)*M1*(pm1-p11.est)+
    (1-M1)*M2*(pm2-p11.est)
  
  Fm1m2
}

trip_nie <- function(dataset,fita,fity,fitm1,fitm2,copula.par)
{
  n <- nrow(dataset)
  #initiate parameters
  A <- dataset$A
  M1 <- dataset$M1
  M2 <- dataset$M2
  Y <- dataset$Y
  
  proba <- fitted(fita)
  da <- dbinom(A,1,proba)
  newdata1 <- newdata0 <- dataset
  newdata1$A <- 1
  newdata0$A <- 0
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  pm1a0 <- expit(predict(fitm1,newdata = newdata0))
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  pm2a0 <- expit(predict(fitm2,newdata = newdata0))
  EYa1 <- predict(fity,newdata = newdata1)
  
  #joint distribution
  fm1 <- med.distrn(pm1a1,pm1a0,pm2a1,pm2a0,1,M1,M2,copula.par = copula.par)
  fm0 <- med.distrn(pm1a1,pm1a0,pm2a1,pm2a0,0,M1,M2,copula.par = copula.par)
  eta10 <- trip_eta.fun(pm1a1,pm1a0,pm2a1,pm2a0,fity,dataset,1,0,copula.par = copula.par)
  eta11 <- trip_eta.fun(pm1a1,pm1a0,pm2a1,pm2a0,fity,dataset,1,1,copula.par = copula.par)
  
  delta0 <- (A/proba)*(fm0/fm1)*(Y-EYa1)+((1-A)/(1-proba))*(EYa1 - eta10) + eta10
  sigma1 <- (A/proba)*(Y - eta11) + eta11
  
  mean(sigma1-delta0)
}

###M1
quad_m1 <- function(dataset,fita,fity,fitm1,fitm2,copula.par)
{
  n <- nrow(dataset)
  #initiate parameters
  A <- dataset$A
  M1 <- dataset$M1
  M2 <- dataset$M2
  Y <- dataset$Y
  dataset$R1 <- 0
  
  proba <- fitted(fita)
  newdata1 <- newdata0 <- dataset
  newdata1$A <- 1
  newdata0$A <- 0
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  pm1a0 <- expit(predict(fitm1,newdata = newdata0))
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  EYa1 <- predict(fity,newdata = newdata1)
  
  one <- delta1.fun(A,Y,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par,EYa1)
  two <- delta2.fun(A,proba,fity,pm1a1,pm1a0,dataset)
  three <- delta3.fun(A,proba,fity,pm2a1,dataset)
  four <- delta4.fun(A,proba,fity,pm1a1,pm1a0,pm2a1,dataset)
  
  delta1 <- mean(one[,1])
  delta2 <- mean(two[,1])
  delta3 <- mean(three[,1])
  delta4 <- mean(four[,1])
  
  delta.quad <- mean(one[,2] + two[,2] + three[,2] + four[,2])
  
  c(delta1,delta2,delta3,delta4,delta.quad)
}

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
  newdata <- dataset
  newdata$A <- 1
  newdata$M1 <- a*pm1a1 + (1-a)*pm1a0
  eta <- predict(fity,newdata = newdata)
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
  newdata <- dataset
  newdata$A <- 1
  newdata$M2 <- pm2a1
  eta <- predict(fity,newdata = newdata)
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
  newdata <- dataset
  newdata$A <- 1
  newdata$M1 <- a*pm1a1 + (1-a)*pm1a0
  newdata$M2 <- pm2a1
  
  gamma <- predict(fity,newdata = newdata)
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



###M2
quad_m2 <- function(dataset,fita,fity,fitm1,fitm2,copula.par)
{
  n <- nrow(dataset)
  #initiate parameters
  A <- dataset$A
  M1 <- dataset$M1
  M2 <- dataset$M2
  Y <- dataset$Y
  dataset$R1 <- 0
  
  proba <- fitted(fita)
  newdata1 <- newdata0 <- dataset
  newdata1$A <- 1
  newdata0$A <- 0
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  pm2a0 <- expit(predict(fitm2,newdata = newdata0))
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  EYa1 <- predict(fity,newdata = newdata1)
  
  one <- delta1.fun2(A,Y,M1,M2,proba,pm2a1,pm2a0,pm1a1,copula.par,EYa1)
  two <- delta2.fun2(A,proba,fity,pm2a1,pm2a0,dataset)
  three <- delta3.fun2(A,proba,fity,pm1a1,dataset)
  four <- delta4.fun2(A,proba,fity,pm2a1,pm2a0,pm1a1,dataset)
  
  delta1 <- mean(one[,1])
  delta2 <- mean(two[,1])
  delta3 <- mean(three[,1])
  delta4 <- mean(four[,1])
  
  delta.quad <- mean(one[,2] + two[,2] + three[,2] + four[,2])
  
  c(delta1,delta2,delta3,delta4,delta.quad)
}

delta1.fun2 <- function(A,Y,M1,M2,proba,pm2a1,pm2a0,pm1a1,copula.par,EYa1)
{
  n <- length(A)
  #joint distribution
  p11.est <- FMC.fun(esttau=copula.par,pc=pm1a1,pm=pm2a1,n)
  Fm1m2 <- M1*M2*p11.est+
    (1-M1)*(1-M2)*(1-pm1a1-pm2a1+p11.est)+
    (1-M2)*M1*(pm1a1-p11.est)+
    (1-M1)*M2*(pm2a1-p11.est)
  
  fm2a1 <- dbinom(M2,1,pm2a1)
  fm2a0 <- dbinom(M2,1,pm2a0)
  fm1a1 <- dbinom(M1,1,pm1a1)
  
  da <- dbinom(A,1,proba)
  ratio <- (A/da)*((fm2a1 - fm2a0)*fm1a1/Fm1m2)
  delta1 <- ratio*Y
  
  comp1 <- ratio*(Y-EYa1)
  
  return(cbind(delta1,comp1))
}

eta1.fun2 <- function(fity,pm2a1,pm2a0,a,dataset)
{
  n <- nrow(dataset)
  newdata <- dataset
  newdata$A <- 1
  newdata$M2 <- a*pm2a1 + (1-a)*pm2a0
  eta <- predict(fity,newdata = newdata)
  eta
}

delta2.fun2 <- function(A,proba,fity,pm2a1,pm2a0,dataset)
{
  n <- nrow(dataset)
  eta.diff <- eta1.fun2(fity,pm2a1,pm2a0,a=1,dataset) - eta1.fun2(fity,pm2a1,pm2a0,a=0,dataset)
  delta2 <- comp2 <- (A/proba)*eta.diff
  return(cbind(delta2,comp2))
}

eta2.fun2 <- function(fity,pm1a1,dataset)
{
  n <- nrow(dataset)
  newdata <- dataset
  newdata$A <- 1
  newdata$M1 <- pm1a1
  eta <- predict(fity,newdata = newdata)
  eta
}

delta3.fun2 <- function(A,proba,fity,pm1a1,dataset)
{
  n <- nrow(dataset)
  da <- dbinom(A,1,proba)
  delta3 <- comp3 <- ((2*A-1)/da)*eta2.fun2(fity,pm1a1,dataset)
  return(cbind(delta3,comp3))
}

gamma1.fun2 <- function(fity,pm2a1,pm2a0,pm1a1,a,dataset)
{
  n <- nrow(dataset)
  newdata <- dataset
  newdata$A <- 1
  newdata$M2 <- a*pm2a1 + (1-a)*pm2a0
  newdata$M1 <- pm1a1
  
  gamma <- predict(fity,newdata = newdata)
  gamma
}

delta4.fun2 <- function(A,proba,fity,pm2a1,pm2a0,pm1a1,dataset)
{
  n <- nrow(dataset)
  da <- dbinom(A,1,proba)
  gamma1 <- gamma1.fun2(fity,pm2a1,pm2a0,pm1a1,1,dataset)
  gamma0 <- gamma1.fun2(fity,pm2a1,pm2a0,pm1a1,0,dataset)
  delta4 <-  gamma1 - gamma0
  comp4 <- delta4 - (2*A/da)*gamma1 + (1/da)*gamma0
  
  return(cbind(delta4,comp4))
}


####stable
ratio.fun1 <- function(A,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par)
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
  
  return(ratio)
}

ratio.fun2 <- function(A,M1,M2,proba,pm2a1,pm2a0,pm1a1,copula.par)
{
  n <- length(A)
  #joint distribution
  p11.est <- FMC.fun(esttau=copula.par,pc=pm1a1,pm=pm2a1,n)
  Fm1m2 <- M1*M2*p11.est+
    (1-M1)*(1-M2)*(1-pm1a1-pm2a1+p11.est)+
    (1-M2)*M1*(pm1a1-p11.est)+
    (1-M1)*M2*(pm2a1-p11.est)
  
  fm2a1 <- dbinom(M2,1,pm2a1)
  fm2a0 <- dbinom(M2,1,pm2a0)
  fm1a1 <- dbinom(M1,1,pm1a1)
  
  da <- dbinom(A,1,proba)
  ratio <- (A/da)*((fm2a1 - fm2a0)*fm1a1/Fm1m2)
  
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
  for.y <- paste0("Y~",paste(c("A","M1","M2","R1",paste0(rep(varname.y,4),1:4)),collapse = "+"))
  
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
  
  R1 <- ratio.fun1(A,M1,M2,proba,pm1a1,pm1a0,pm2a1,copula.par=esttau)
  dataset$R1 <- as.numeric(R1)
  
  fity <- lm(for.y,data = dataset,weights=weighta)
  quad.stable <- quad_m1(dataset,fita,fity,fitm1,fitm2,esttau)
  
  return(quad.stable[5])
}

stable1.m2 <- function(dataset,fita,M1T=T,M2T=T,YT=T)
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
  for.y <- paste0("Y~",paste(c("A","M1","M2","R1",paste0(rep(varname.y,4),1:4)),collapse = "+"))
  
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
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  pm2a0 <- expit(predict(fitm2,newdata = newdata0))
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  
  R1 <- ratio.fun2(A,M1,M2,proba,pm2a1,pm2a0,pm1a1,copula.par=esttau)
  dataset$R1 <- as.numeric(R1)
  
  fity <- lm(for.y,data = dataset,weights=weighta)
  quad.stable <- quad_m2(dataset,fita,fity,fitm1,fitm2,esttau)
  
  return(quad.stable[5])
} 

stable1.nie <- function(dataset,fita,M1T=T,M2T=T,YT=T)
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
  for.y <- paste0("Y~",paste(c("A","M1","M2",paste0(rep(varname.y,4),1:4)),collapse = "+"))
  
  proba <- fitted(fita)
  dataset$weighta <- (A*proba+(1-A)*(1-proba))^(-1)
  fitm1 <- glm(for.m1,family = binomial(link = "logit"),data = dataset,
               weights = weighta)
  fitm2 <- glm(for.m2,family = binomial(link = "logit"),data = dataset,
               weights = weighta)
  
  copula.par <- estp11(fitm2,fitm1,M2,M1,nrow(dataset),dataset)
  newdata1 <- newdata0 <- dataset
  newdata1$A <- 1
  newdata0$A <- 0
  pm1a1 <- expit(predict(fitm1,newdata = newdata1))
  pm1a0 <- expit(predict(fitm1,newdata = newdata0))
  pm2a1 <- expit(predict(fitm2,newdata = newdata1))
  pm2a0 <- expit(predict(fitm2,newdata = newdata0))
  
  #joint distribution
  fm1 <- med.distrn(pm1a1,pm1a0,pm2a1,pm2a0,1,M1,M2,copula.par = copula.par)
  fm0 <- med.distrn(pm1a1,pm1a0,pm2a1,pm2a0,0,M1,M2,copula.par = copula.par)
  
  dataset$weighty <- fm0/(proba*fm1)
  fity <- lm(for.y,data = dataset,weights=weighty)
  trip.stable.nie <- trip_nie(dataset,fita,fity,fitm1,fitm2,copula.par)
  
  trip.stable.nie
}


##########main

#helper functions
N <- 5000
#parameters
tau_0 <- 1.2
tau_A <- 0
tau_X <- c(0,0,0,0)
#baseline covariates
#parameters for A
aU <- c(-1,0.5,-0.25,-0.1) 
#marginal for mediators
#M1 marginal parameters
m10 <- -1.6
m1A <- 2
m1U <- c(1,-0.8,0.6,-1)


#M2 marginal parameters
m20 <- -1.5
m2A <- 2
m2U <- c(1,-0.5,0.9,-1)


#generate Y
y0 <- 210
yA <- 1
yM2 <- 1
yM1 <- -50
yU <- c(27.4,13.7,13.7,13.7)
sdy <- 50

#calculating the true value exactly
U1 <- U2 <- U3 <- U4 <- rep(NA,N)
for(i in 1:N)
{
  U.temp <- mvrnorm(1,mu=cbind(0,0,0,0), Sigma=diag(c(1,1,1,1)))
  U1[i] <- U.temp[1]
  U2[i] <- U.temp[2]
  U3[i] <- U.temp[3]
  U4[i] <- U.temp[4]
}
U.row <- cbind(U1,U2,U3,U4)
U <- t(U.row)


Ey <- function(a,m2,m1)
{
  return(y0+yA*a+yM2*m2+yM1*m1+yU %*% U)
}

gamma1_1 <- mean(Ey(1,1,1)*expit(m20+m2A+m2U%*%U)*expit(m10+m1A+m1U%*%U)+
                   Ey(1,1,0)*expit(m20+m2A+m2U%*%U)*(1-expit(m10+m1A+m1U%*%U))+
                   Ey(1,0,1)*(1-expit(m20+m2A+m2U%*%U))*expit(m10+m1A+m1U%*%U)+
                   Ey(1,0,0)*(1-expit(m20+m2A+m2U%*%U))*(1-expit(m10+m1A+m1U%*%U))
)

gamma1_0 <- mean(Ey(1,1,1)*expit(m20+m2A+m2U%*%U)*expit(m10+m1U%*%U)+
                   Ey(1,1,0)*expit(m20+m2A+m2U%*%U)*(1-expit(m10+m1U%*%U))+
                   Ey(1,0,1)*(1-expit(m20+m2A+m2U%*%U))*expit(m10+m1U%*%U)+
                   Ey(1,0,0)*(1-expit(m20+m2A+m2U%*%U))*(1-expit(m10+m1U%*%U))
)

deltaM1 <- gamma1_1 - gamma1_0

gamma2_1 <- mean(Ey(1,1,1)*expit(m10+m1A+m1U%*%U)*expit(m20+m2A+m2U%*%U)+
                   Ey(1,0,1)*expit(m10+m1A+m1U%*%U)*(1-expit(m20+m2A+m2U%*%U))+
                   Ey(1,1,0)*(1-expit(m10+m1A+m1U%*%U))*expit(m20+m2A+m2U%*%U)+
                   Ey(1,0,0)*(1-expit(m10+m1A+m1U%*%U))*(1-expit(m20+m2A+m2U%*%U))
)

gamma2_0 <- mean(Ey(1,1,1)*expit(m10+m1A+m1U%*%U)*expit(m20+m2U%*%U)+
                   Ey(1,0,1)*expit(m10+m1A+m1U%*%U)*(1-expit(m20+m2U%*%U))+
                   Ey(1,1,0)*(1-expit(m10+m1A+m1U%*%U))*expit(m20+m2U%*%U)+
                   Ey(1,0,0)*(1-expit(m10+m1A+m1U%*%U))*(1-expit(m20+m2U%*%U))
)
deltaM2 <- gamma2_1 - gamma2_0

p11.true_a1 <- FMC.fun(esttau = tau_0,pc=expit(m10+m1A+m1U%*%U),pm=expit(m20+m2A+m2U%*%U),N)
p10.true_a1 <- expit(m10+m1A+m1U%*%U) - p11.true_a1
p01.true_a1 <- expit(m20+m2A+m2U%*%U) - p11.true_a1
p00.true_a1 <- 1 - expit(m10+m1A+m1U%*%U) - expit(m20+m2A+m2U%*%U) + p11.true_a1

rho1 <- mean(Ey(1,1,1)*p11.true_a1+ Ey(1,0,1)*p10.true_a1 + Ey(1,1,0)*p01.true_a1 +
               Ey(1,0,0)*p00.true_a1)

p11.true_a0 <- FMC.fun(esttau = tau_0,pc=expit(m10+m1U%*%U),pm=expit(m20+m2U%*%U),N)
p10.true_a0 <- expit(m10+m1U%*%U) - p11.true_a0
p01.true_a0 <- expit(m20+m2U%*%U) - p11.true_a0
p00.true_a0 <- 1 - expit(m10+m1U%*%U) - expit(m20+m2U%*%U) + p11.true_a0

rho0 <- mean(Ey(1,1,1)*p11.true_a0 + Ey(1,0,1)*p10.true_a0 + Ey(1,1,0)*p01.true_a0 +
               Ey(1,0,0)*p00.true_a0)

nie.true <- rho1 - rho0

int.true <- nie.true - deltaM1 - deltaM2


######use
bootstrap <- function(function.use,dataset)
{
  
  nsize <- nrow(dataset)
  bsample.index <- sample(nsize,nsize,replace = T)
  
  bdata <- dataset[bsample.index,]
  res.one <- function.use(bdata)
  res.one
}
corr.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+M1+M2,data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- stable1.m1(dataset.Z,fita,M1T=T,M2T=T,YT=T)
  
  #######EIEM2#
  res.m2 <- quad_m2(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m2 <- stable1.m2(dataset.Z,fita,M1T=T,M2T=T,YT=T)
  
  #########nie
  res.nie <- trip_nie(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.nie <- stable1.nie(dataset.Z,fita,M1T=T,M2T=T,YT=T)
  
  #########int
  res.int <- res.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  stable.int <- stable.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  c(res.m1,stable.m1,res.m2,stable.m2,res.int,stable.int)
}

prop.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+M1+M2,data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- stable1.m1(dataset.Z,fita,M1T=T,M2T=T,YT=T)
  
  #######EIEM2#
  res.m2 <- quad_m2(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m2 <- stable1.m2(dataset.Z,fita,M1T=T,M2T=T,YT=T)
  
  #########nie
  res.nie <- trip_nie(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.nie <- stable1.nie(dataset.Z,fita,M1T=T,M2T=T,YT=T)
  
  #########int
  res.int <- res.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  stable.int <- stable.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  c(res.m1,stable.m1,res.m2,stable.m2,res.int,stable.int)
}
outcome.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+M1+M2,data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- stable1.m1(dataset.Z,fita,M1T=T,M2T=T,YT=F)
  
  #######EIEM2#
  res.m2 <- quad_m2(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m2 <- stable1.m2(dataset.Z,fita,M1T=T,M2T=T,YT=F)
  
  #########nie
  res.nie <- trip_nie(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.nie <- stable1.nie(dataset.Z,fita,M1T=T,M2T=T,YT=F)
  
  #########int
  res.int <- res.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  stable.int <- stable.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  c(res.m1,stable.m1,res.m2,stable.m2,res.int,stable.int)
}

medone.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+M1+M2,data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- stable1.m1(dataset.Z,fita,M1T=F,M2T=T,YT=T)
  
  #######EIEM2#
  res.m2 <- quad_m2(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m2 <- stable1.m2(dataset.Z,fita,M1T=F,M2T=T,YT=T)
  
  #########nie
  res.nie <- trip_nie(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.nie <- stable1.nie(dataset.Z,fita,M1T=F,M2T=T,YT=T)
  
  #########int
  res.int <- res.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  stable.int <- stable.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  c(res.m1,stable.m1,res.m2,stable.m2,res.int,stable.int)
}

medtwo.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+U1+U2+U3+U4+M1+M2,data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- stable1.m1(dataset.Z,fita,M1T=T,M2T=F,YT=T)
  
  #######EIEM2#
  res.m2 <- quad_m2(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m2 <- stable1.m2(dataset.Z,fita,M1T=T,M2T=F,YT=T)
  
  #########nie
  res.nie <- trip_nie(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.nie <- stable1.nie(dataset.Z,fita,M1T=T,M2T=F,YT=T)
  
  #########int
  res.int <- res.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  stable.int <- stable.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  c(res.m1,stable.m1,res.m2,stable.m2,res.int,stable.int)
}
wrong.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- lm(Y~A+Z1+Z2+Z3+Z4+M1+M2,data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- stable1.m1(dataset.Z,fita,M1T=F,M2T=F,YT=F)
  
  #######EIEM2#
  res.m2 <- quad_m2(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m2 <- stable1.m2(dataset.Z,fita,M1T=F,M2T=F,YT=F)
  
  #########nie
  res.nie <- trip_nie(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.nie <- stable1.nie(dataset.Z,fita,M1T=F,M2T=F,YT=F)
  
  #########int
  res.int <- res.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  stable.int <- stable.nie - c(res.m1,stable.m1) - c(res.m2[1],res.m2[3],res.m2[2],res.m2[4],res.m2[5],stable.m2)
  c(res.m1,stable.m1,res.m2,stable.m2,res.int,stable.int)
}


