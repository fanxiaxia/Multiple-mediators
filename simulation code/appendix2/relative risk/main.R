##############
##########main

#set parameters and calculate "true" value
N <- 20000
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
y0 <- 1
yA <- -1.6
yM2 <- -1
yM1 <- -1.2
yAM1 <- -0.2
yAM2 <- -0.1
yM1M2 <- 0.05
yU <- c(1,-0.2,-0.8,-1)



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
  return(exp(y0+yA*a+yM2*m2+yM1*m1+yAM1*a*m1+yAM2*a*m2 + yM1M2*m1*m2 + yU %*% U))
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

######use the following functions to compute either
#1. moment estimators and original quadruply robust estimator
#2. the stabilized quadruply robust estimator
######for different incorrectly specified nuisance models.

corr.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- glm(Y~A+U1+U2+U3+U4+M1+M2+A*M1+A*M2+M1*M2,family = poisson(link = "log"),data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- tryCatch(stable1.m1(dataset.Z,fita,M1T=T,M2T=T,YT=T,esttau),error = function(e){print("bad dataset!");rep(NA,2)})
  c(res.m1,stable.m1)
}

prop.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- glm(Y~A+U1+U2+U3+U4+M1+M2+A*M1+A*M2+M1*M2,family = poisson(link = "log"),data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- tryCatch(stable1.m1(dataset.Z,fita,M1T=T,M2T=T,YT=T,esttau),error = function(e){print("bad dataset!");rep(NA,2)})
  c(res.m1,stable.m1)
}

outcome.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- glm(Y~A+Z1+Z2+Z3+Z4+M1+M2+A*M1+A*M2+M1*M2,family = poisson(link = "log"),data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- tryCatch(stable1.m1(dataset.Z,fita,M1T=T,M2T=T,YT=F,esttau),error = function(e){print("bad dataset!");rep(NA,2)})
  c(res.m1,stable.m1)
}

medone.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- glm(Y~A+U1+U2+U3+U4+M1+M2+A*M1+A*M2+M1*M2,family = poisson(link = "log"),data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- tryCatch(stable1.m1(dataset.Z,fita,M1T=F,M2T=T,YT=T,esttau),error = function(e){print("bad dataset!");rep(NA,2)})
  c(res.m1,stable.m1)
}

medtwo.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~U1+U2+U3+U4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+U1+U2+U3+U4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- glm(Y~A+U1+U2+U3+U4+M1+M2+A*M1+A*M2+M1*M2,family = poisson(link = "log"),data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- tryCatch(stable1.m1(dataset.Z,fita,M1T=T,M2T=F,YT=T,esttau),error = function(e){print("bad dataset!");rep(NA,2)})
  c(res.m1,stable.m1)
  
}
wrong.fun <- function(dataset.Z)
{
  n <- nrow(dataset.Z)
  fita <- glm(A~Z1+Z2+Z3+Z4-1,family = binomial(link = "logit"),data = dataset.Z)
  fitm2.mar <- glm(M2~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fitm1.mar <- glm(M1~A+Z1+Z2+Z3+Z4,family = binomial(link = "logit"),data = dataset.Z)
  fity <- glm(Y~A+Z1+Z2+Z3+Z4+M1+M2+A*M1+A*M2+M1*M2,family = poisson(link = "log"),data = dataset.Z)
  
  esttau <- estp11(fitm2.mar,fitm1.mar,dataset.Z$M2,dataset.Z$M1,n,dataset.Z)
  
  #######EIEM1#  
  res.m1 <- quad_m1(dataset=dataset.Z,fita,fity,fitm1=fitm1.mar,fitm2=fitm2.mar,copula.par=esttau)
  stable.m1 <- tryCatch(stable1.m1(dataset.Z,fita,M1T=F,M2T=F,YT=F,esttau),error = function(e){print("bad dataset!");rep(NA,2)})
  c(res.m1,stable.m1)
}

###replicate
truevector <- c(rep(deltaM1,6))
do.one <- function(size,type)
{
  #  set.seed(index+233)
  dataset.Z <- data_gen(size)
  function.use <- match.fun(paste0(type,".fun"))
  res <- function.use(dataset.Z)
  res
}

verify <- function(i){
  
  res <- do.one(size,type)
  
  return(res)
}

#the stabilized algorithm may not converge, so B may need to be set to larger number for stabilized estimator.
B <- 1000
trials <- 1:B

size <- 1000
type.vector <- c("corr","prop","medone","medtwo","outcome","wrong")
for(type in type.vector)
{
  res.temp <- mclapply(trials,verify,mc.cores=detectCores())
}


