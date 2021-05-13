setwd("/Users/fanxia/Dropbox/kang-new simulation-biometrika/binarymediator/final estimation/")
source("bpsource.R")
library(parallel)

truevector <- c(rep(deltaM1,6),rep(deltaM2,6),rep(int.true,12))
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

trials <- 1:1000

size.vector <- c(2500,5000)
type.vector <- c("corr","prop","medone","medtwo","outcome","wrong")

for(size in size.vector)
{
  biastable <- sdtable <- matrix(NA,nrow = length(type.vector),ncol = 6*4)
  rownames(biastable) <- rownames(sdtable) <- type.vector
  colnames(biastable) <- colnames(sdtable) <- c(rep("deltaM1",6),rep("deltaM2",6),rep("int",12))
  for(type in type.vector)
  {
    res.temp <- mclapply(trials,verify,mc.cores=5)
    res.matrix.temp <- matrix(unlist(res.temp),nrow=24)
    biastable[type,] <- apply(res.matrix.temp,1,mean) - truevector
    sdtable[type,] <- apply(res.matrix.temp,1,sd)
  }
  write.csv(biastable,paste0(size,"bias.csv"))
  write.csv(sdtable,paste0(size,"sd.csv"))
}

