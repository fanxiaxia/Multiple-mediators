library(snowfall)

do.one.np <- function(size,i) 
{
  set.seed(i+233)
  dataset.Z <- data_gen(size)
  res <- npest(dataset.Z,V=2)
  res
}


verify <- function(i){
  
  res <- do.one.np(size,i)
  
  return(res)
}

set.seed(9012)
M <- 1000

size_seq <- c(250,500,1000,2500,5000)

for(size in size_seq)
{
  sfExport("size") 
  res <- sfSapply(1:M, verify)
}