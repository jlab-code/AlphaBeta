


num.Core <- function(NumFiles){
  # no_cores <- detectCores()
  # if (no_cores > NumFiles){
  #   no_cores <- NumFiles
  # }else{
  #   no_cores <- no_cores - 1
  # }
  cl <- makeCluster(NumFiles, type="FORK")
  registerDoParallel(cl)
  return(cl)
}

# take 5 digit of decimal value posteriorMax column
floorDec<-function(rc.Mean ,x){
  y <- function(x, level=1) round(x - 5*10^(-level-1), level)
  res <-y(as.numeric(rc.Mean),x)
  return(res)
}
