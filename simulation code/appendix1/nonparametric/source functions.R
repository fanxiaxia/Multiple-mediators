library(SuperLearner)
library(caret)
library(MASS)
###helper
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
FMC.fun <- function(esttau,pc,pm)
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

###############true value
set.seed(336)
#setwd("/Users/fanxia/Desktop/new simulation/new binary MC/")
#helper functions
n <- 10000
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
U1 <- U2 <- U3 <- U4 <- rep(NA,n)
for(i in 1:n)
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

p11.true_a1 <- FMC.fun(esttau = tau_0,pc=expit(m10+m1A+m1U%*%U),pm=expit(m20+m2A+m2U%*%U))
p10.true_a1 <- expit(m10+m1A+m1U%*%U) - p11.true_a1
p01.true_a1 <- expit(m20+m2A+m2U%*%U) - p11.true_a1
p00.true_a1 <- 1 - expit(m10+m1A+m1U%*%U) - expit(m20+m2A+m2U%*%U) + p11.true_a1

rho1 <- mean(Ey(1,1,1)*p11.true_a1+ Ey(1,0,1)*p10.true_a1 + Ey(1,1,0)*p01.true_a1 +
               Ey(1,0,0)*p00.true_a1)

p11.true_a0 <- FMC.fun(esttau = tau_0,pc=expit(m10+m1U%*%U),pm=expit(m20+m2U%*%U))
p10.true_a0 <- expit(m10+m1U%*%U) - p11.true_a0
p01.true_a0 <- expit(m20+m2U%*%U) - p11.true_a0
p00.true_a0 <- 1 - expit(m10+m1U%*%U) - expit(m20+m2U%*%U) + p11.true_a0

rho0 <- mean(Ey(1,1,1)*p11.true_a0 + Ey(1,0,1)*p10.true_a0 + Ey(1,1,0)*p01.true_a0 +
               Ey(1,0,0)*p00.true_a0)

nie.true <- rho1 - rho0
int.true <- nie.true - deltaM1 - deltaM2
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
###########
npest<- function(dataset,V=2)
{
  index <- createFolds(1:nrow(dataset),V)
  int <- nie <- quadm1 <- quadm2 <- NULL
  for(v in 1:V)
  {
    trainingset <- dataset[-index[[v]],]
    predictset <- dataset[index[[v]],]
    
    sl_a <- SuperLearner(Y = trainingset$A, X = trainingset[,c("Z1","Z2","Z3","Z4")],
                         family = binomial(), SL.library = c("SL.glm", "SL.randomForest"),
                         cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL))
    sl_a_con <- SuperLearner(Y = trainingset$A, X = trainingset[,c("Z1","Z2","Z3","Z4","M1","M2")],
                             family = binomial(), SL.library = c("SL.glm", "SL.randomForest"),
                             cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL))
    sl_m2_mar <- SuperLearner(Y = trainingset$M2, X = trainingset[,c("A","Z1","Z2","Z3","Z4")],
                              family = binomial(), SL.library = c("SL.glm", "SL.randomForest"),
                              cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL))
    sl_m1_mar <- SuperLearner(Y = trainingset$M1, X = trainingset[,c("A","Z1","Z2","Z3","Z4")],
                              family = binomial(), SL.library = c("SL.glm", "SL.randomForest"),
                              cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL))
    sl_m2_con <- SuperLearner(Y = trainingset$M2, X = trainingset[,c("A","M1","Z1","Z2","Z3","Z4")],
                              family = binomial(), SL.library = c("SL.glm", "SL.randomForest"),
                              cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL))
    sl_m1_con <- SuperLearner(Y = trainingset$M1, X = trainingset[,c("A","M2","Z1","Z2","Z3","Z4")],
                              family = binomial(), SL.library = c("SL.glm", "SL.randomForest"),
                              cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL))
    sl_y <- SuperLearner(Y = trainingset$Y, X = trainingset[,c("A", "M1","M2", "Z1","Z2","Z3","Z4")],
                         family     = gaussian(),SL.library = c("SL.glm", "SL.randomForest"),
                         cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL))
    
    traina1 <- trainingset
    traina1$A <- 1
    Eya1.train <- predict(sl_y,newdata = traina1,onlySL = T)$pred
    
    sl_Ey1 <- SuperLearner(Y = Eya1.train, X = trainingset[,c("A","Z1","Z2","Z3","Z4")],
                           family     = gaussian(),SL.library = c("SL.glm", "SL.randomForest"),
                           cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL))
    
    nietrip <- function(sl_a,sl_y,sl_Ey1,sl_a_con,predictset)
    {
      proba <- predict(sl_a, newdata = predictset, onlySL = T)$pred
      datam1 <- datam0 <-  predictset
      datam1$A <- 1  
      datam0$A <- 0
      EYa1 <- predict(sl_y, newdata = datam1, onlySL = TRUE)$pred
      ################    NIE    #############    
      eta10 <- predict(sl_Ey1,newdata = datam0, onlySL = T)$pred
      eta11 <- predict(sl_Ey1,newdata = datam1, onlySL = T)$pred
      
      proba_con <- predict(sl_a_con, newdata = predictset, onlySL = T)$pred
      ratio.nie <- (1 - proba_con)/proba_con
      
      delta0.1 <- (predictset$A/(1-proba))*(ratio.nie)*(predictset$Y-EYa1)+((1-predictset$A)/(1 -proba))*(EYa1 - eta10) + eta10
      sigma1.1 <- (predictset$A/proba)*(predictset$Y - eta11) + eta11
      
      (sigma1.1 - delta0.1)
    }
    nie.temp <- nietrip(sl_a,sl_y,sl_Ey1,sl_a_con,predictset)
    nie <- c(nie,nie.temp)
    ########################################
    ################    M1    ##############      
    m1quad <- function(sl_a,sl_m1_mar,sl_m1_con,sl_m2_mar,sl_y,predictset)
    {
      proba <- dbinom(predictset$A,1,predict(sl_a, newdata = predictset, onlySL = T)$pred)
      datam1 <- datam0 <-  predictset
      datam1$A <- 1  
      datam0$A <- 0
      
      probm1_1 <- predict(sl_m1_mar, newdata = datam1, onlySL = T)$pred
      probm1_0 <- predict(sl_m1_mar, newdata = datam0, onlySL = T)$pred
      probm1_1_cond <- predict(sl_m1_con, newdata = datam1, onlySL = T)$pred
      probm2_1 <- predict(sl_m2_mar, newdata = datam1, onlySL = T)$pred
      
      dm10x <- dbinom(predictset$M1,1,probm1_0)
      dm11x <- dbinom(predictset$M1,1,probm1_1)
      dm11x.con <- dbinom(predictset$M1,1,probm1_1_cond)
      
      datay1 <- predictset 
      datay1$A <- 1
      EY1 = predict(sl_y, newdata = datay1, onlySL = TRUE)$pred
      
      ratio1 <- (predictset$A/proba)*(dm11x-dm10x)/dm11x.con
      est_comp1 <- ratio1*(predictset$Y-EY1)
      
      etafun_1 <- function(pm10,pm11)
      {
        dataeta1_1 <- dataeta1_0 <- predictset
        dataeta1_1$A <- dataeta1_0$A <- 1
        dataeta1_1$M1 <- 1
        dataeta1_0$M1 <- 0
        eta_1 <- predict(sl_y, newdata = dataeta1_1, onlySL = TRUE)$pred*pm11 + 
          predict(sl_y, newdata = dataeta1_0, onlySL = TRUE)$pred*(1 - pm11)
        
        eta_0 <- predict(sl_y, newdata = dataeta1_1, onlySL = TRUE)$pred*pm10 + 
          predict(sl_y, newdata = dataeta1_0, onlySL = TRUE)$pred*(1 - pm10)
        
        eta <- eta_1 - eta_0
      }
      
      ratio2 <- predictset$A/proba
      est_delta2 <- est_comp2 <- ratio2*etafun_1(pm11 = probm1_1, pm10 = probm1_0)
      
      etafun_2 <-function(pm2)
      {
        dataeta2_1 <- dataeta2_0 <- predictset
        dataeta2_1$A <- 1
        dataeta2_1$M2 <- 1
        dataeta2_0$M2 <- 0
        eta <- predict(sl_y, newdata = dataeta2_1, onlySL = TRUE)$pred*pm2 + 
          predict(sl_y, newdata = dataeta2_0, onlySL = TRUE)$pred*(1-pm2)
      }
      
      ratio3 <- (2*predictset$A-1)/proba
      est_delta3 <- est_comp3 <- ratio3*etafun_2(pm2=probm2_1)
      
      gammafun <- function(X,pm10,pm11,pm2)
      {
        datac1m1 <- datac1m0 <- datac0m0 <- datac0m1 <- predictset
        datac1m1$A <- datac0m0$A <- datac1m0$A <- datac0m1$A <- 1
        
        datac1m1$M1 <- datac1m1$M2 <- 1
        datac1m0$M1 <- datac0m1$M2 <- 1
        datac1m0$M2 <- datac0m1$M1 <- 0
        
        
        tau_1 <- predict(sl_y, newdata = datac1m1, onlySL = TRUE)$pred*pm11*pm2 +
          predict(sl_y, newdata = datac0m1, onlySL = TRUE)$pred*(1-pm11)*pm2 +
          predict(sl_y, newdata = datac1m0, onlySL = TRUE)$pred*pm11*(1-pm2) +
          predict(sl_y, newdata = datac0m0, onlySL = TRUE)$pred*(1-pm11)*(1-pm2)
        
        tau_0 <- predict(sl_y, newdata = datac1m1, onlySL = TRUE)$pred*pm10*pm2 +
          predict(sl_y, newdata = datac0m1, onlySL = TRUE)$pred*(1-pm10)*pm2 +
          predict(sl_y, newdata = datac1m0, onlySL = TRUE)$pred*pm10*(1-pm2) +
          predict(sl_y, newdata = datac0m0, onlySL = TRUE)$pred*(1-pm10)*(1-pm2)
        
        tau <- tau_1 - tau_0
      }
      est_delta4 <- gammafun(pm10=probm1_0,pm11=probm1_1,pm2=probm2_1)
      est_comp4 <- est_delta4
      
      iden.term <- function(a,pm1,pc0,pc1)
      {
        pM1 <- pm1
        pCa <- a*pc1 + (1 - a)*pc0
        
        datac1m1 <- datac1m0 <- datac0m0 <- datac0m1 <- predictset
        datac1m1$A <- datac0m0$A <- datac1m0$A <- datac0m1$A <- 1
        
        datac1m1$M1 <- datac1m1$M2 <- 1
        datac1m0$M1 <- datac0m1$M2 <- 1
        datac1m0$M2 <- datac0m1$M1 <- 0
        
        
        thi <- predict(sl_y, newdata = datac1m1, onlySL = TRUE)$pred*pM1*pCa+
          predict(sl_y, newdata = datac1m0, onlySL = TRUE)$pred*(1-pM1)*pCa+
          predict(sl_y, newdata = datac0m1, onlySL = TRUE)$pred*(1-pCa)*pM1+
          predict(sl_y, newdata = datac0m0, onlySL = TRUE)$pred*(1-pCa)*(1-pM1)
        
        thi
      }
      
      est_comp5 <- (-2*predictset$A/proba)*iden.term(1,pm1=probm2_1,pc0=probm1_0,pc1=probm1_1)-
        (-predictset$A/proba-(1-predictset$A)/proba)*iden.term(0,pm1=probm2_1,pc0=probm1_0,pc1=probm1_1)
      
      est_comp1 + est_comp2 + est_comp3 + est_comp4 + est_comp5
    }
    m1.temp <- m1quad(sl_a,sl_m1_mar,sl_m1_con,sl_m2_mar,sl_y,predictset)
    quadm1 <- c(quadm1, m1.temp)
    #############################################
    #####################  M2  ##################
    m2quad <- function(sl_a,sl_m2_mar,sl_m2_con,sl_m1_mar,sl_y,predictset)
    {
      proba <- dbinom(predictset$A,1,predict(sl_a, newdata = predictset, onlySL = T)$pred)
      datam1 <- datam0 <-  predictset
      datam1$A <- 1  
      datam0$A <- 0
      
      probm2_1 <- predict(sl_m2_mar, newdata = datam1, onlySL = T)$pred
      probm2_0 <- predict(sl_m2_mar, newdata = datam0, onlySL = T)$pred
      probm2_1_cond <- predict(sl_m2_con, newdata = datam1, onlySL = T)$pred
      probm1_1 <- predict(sl_m1_mar, newdata = datam1, onlySL = T)$pred
      
      dm20x <- dbinom(predictset$M2,1,probm2_0)
      dm21x <- dbinom(predictset$M2,1,probm2_1)
      dm21x.con <- dbinom(predictset$M2,1,probm2_1_cond)
      
      datay1 <- predictset 
      datay1$A <- 1
      EY1 = predict(sl_y, newdata = datay1, onlySL = TRUE)$pred
      
      ratio1 <- (predictset$A/proba)*(dm21x-dm20x)/dm21x.con
      est_comp1 <- ratio1*(predictset$Y-EY1)
      
      etafun_1 <- function(pm20,pm21)
      {
        dataeta1_1 <- dataeta1_0 <- predictset
        dataeta1_1$A <- dataeta1_0$A <- 1
        dataeta1_1$M2 <- 1
        dataeta1_0$M2 <- 0
        eta_1 <- predict(sl_y, newdata = dataeta1_1, onlySL = TRUE)$pred*pm21 + 
          predict(sl_y, newdata = dataeta1_0, onlySL = TRUE)$pred*(1 - pm21)
        
        eta_0 <- predict(sl_y, newdata = dataeta1_1, onlySL = TRUE)$pred*pm20 + 
          predict(sl_y, newdata = dataeta1_0, onlySL = TRUE)$pred*(1 - pm20)
        
        eta <- eta_1 - eta_0
      }
      
      ratio2 <- predictset$A/proba
      est_delta2 <- est_comp2 <- ratio2*etafun_1(pm21 = probm2_1, pm20 = probm2_0)
      
      etafun_2 <-function(pm1)
      {
        dataeta2_1 <- dataeta2_0 <- predictset
        dataeta2_1$A <- 1
        dataeta2_1$M1 <- 1
        dataeta2_0$M1 <- 0
        eta <- predict(sl_y, newdata = dataeta2_1, onlySL = TRUE)$pred*pm1 + 
          predict(sl_y, newdata = dataeta2_0, onlySL = TRUE)$pred*(1-pm1)
      }
      
      ratio3 <- (2*predictset$A-1)/proba
      est_delta3 <- est_comp3 <- ratio3*etafun_2(pm1=probm1_1)
      
      gammafun <- function(X,pm20,pm21,pm1)
      {
        datac1m1 <- datac1m0 <- datac0m0 <- datac0m1 <- predictset
        datac1m1$A <- datac0m0$A <- datac1m0$A <- datac0m1$A <- 1
        
        datac1m1$M1 <- datac1m1$M2 <- 1
        datac1m0$M2 <- datac0m1$M1 <- 1
        datac1m0$M1 <- datac0m1$M2 <- 0
        
        
        tau_1 <- predict(sl_y, newdata = datac1m1, onlySL = TRUE)$pred*pm21*pm1 +
          predict(sl_y, newdata = datac0m1, onlySL = TRUE)$pred*(1-pm21)*pm1 +
          predict(sl_y, newdata = datac1m0, onlySL = TRUE)$pred*pm21*(1-pm1) +
          predict(sl_y, newdata = datac0m0, onlySL = TRUE)$pred*(1-pm21)*(1-pm1)
        
        tau_0 <- predict(sl_y, newdata = datac1m1, onlySL = TRUE)$pred*pm20*pm1 +
          predict(sl_y, newdata = datac0m1, onlySL = TRUE)$pred*(1-pm20)*pm1 +
          predict(sl_y, newdata = datac1m0, onlySL = TRUE)$pred*pm20*(1-pm1) +
          predict(sl_y, newdata = datac0m0, onlySL = TRUE)$pred*(1-pm20)*(1-pm1)
        
        tau <- tau_1 - tau_0
      }
      est_delta4 <- gammafun(pm20=probm2_0,pm21=probm2_1,pm1=probm1_1)
      est_comp4 <- est_delta4
      
      iden.term <- function(a,pm1,pc0,pc1)
      {
        pM1 <- pm1
        pCa <- a*pc1 + (1 - a)*pc0
        
        datac1m1 <- datac1m0 <- datac0m0 <- datac0m1 <- predictset
        datac1m1$A <- datac0m0$A <- datac1m0$A <- datac0m1$A <- 1
        
        datac1m1$M1 <- datac1m1$M2 <- 1
        datac1m0$M2 <- datac0m1$M1 <- 1
        datac1m0$M1 <- datac0m1$M2 <- 0
        
        
        thi <- predict(sl_y, newdata = datac1m1, onlySL = TRUE)$pred*pM1*pCa+
          predict(sl_y, newdata = datac1m0, onlySL = TRUE)$pred*(1-pM1)*pCa+
          predict(sl_y, newdata = datac0m1, onlySL = TRUE)$pred*(1-pCa)*pM1+
          predict(sl_y, newdata = datac0m0, onlySL = TRUE)$pred*(1-pCa)*(1-pM1)
        
        thi
      }
      
      est_comp5 <- (-2*predictset$A/proba)*iden.term(1,pm1=probm1_1,pc0=probm2_0,pc1=probm2_1)-
        (-predictset$A/proba-(1-predictset$A)/proba)*iden.term(0,pm1=probm1_1,pc0=probm2_0,pc1=probm2_1)
      
      
      (est_comp1 + est_comp2 + est_comp3 + est_comp4 + est_comp5)
    }
    m2.temp <- m2quad(sl_a,sl_m2_mar,sl_m2_con,sl_m1_mar,sl_y,predictset)
    quadm2 <- c(quadm2,m2.temp)
    #############################################  
    ############### INT #########################
    int.temp <- nie.temp - m1.temp - m2.temp
    int <- c(int,int.temp)
  }      
  est.m1 <- mean(quadm1)
  sd.m1 <- sd(quadm1)
  est.m2 <- mean(quadm2)
  sd.m2 <- sd(quadm2)
  est.nie <- mean(nie)
  sd.nie <- sd(nie)
  est.int <- mean(int)
  sd.int <- sd(int)
  
  return(c(est.m1,sd.m1,est.m2,sd.m2,est.nie,sd.nie,est.int,sd.int))
}