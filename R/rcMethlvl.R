#' Calculating rc.Meth.lvl
#'
#' Estimating epimutation rates from high-throughput DNA methylation data
#' @param genTable Generation table name, you can find sample file in
#' "extdata" called "generations.fn"
#' @param cytosine Type of cytosine (CHH/CHG/CG)
#' @param posteriorMaxFilter Filter value, based on posteriorMax
#' @param nThread number of threads, default is 2.
#' @import      dplyr
#' @import      parallel
#' @import      doParallel
#' @import      foreach
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom  data.table as.data.table
#' @importFrom  stringr str_replace_all
#' @importFrom  foreach %dopar%
#' @importFrom  doParallel registerDoParallel
#' @importFrom  gtools mixedorder
#' @return rc meth lvl.
#' @export
#' @examples
#'## Get some toy data
#' file <- system.file("extdata","tm_generations.fn", package="AlphaBeta")
#' rc.meth.lvl(file, "CG", 0.99)


# running join
rc.meth.lvl <- function(genTable, cytosine, posteriorMaxFilter, nThread=2){

  inputCheck(genTable, cytosine, posteriorMaxFilter)
  genTable <- fread(genTable)
  mt <- startTime("Preparing data-sets...\n")
  i<-NULL

  if(.Platform$OS.type == "unix") {
  cal<-num.Core(nThread)
  list.rc <- foreach(i=seq_len(length(genTable$filename))) %dopar% rcRun(genTable$filename[i],cytosine, posteriorMaxFilter, genTable)
  stopCluster(cal)
  }else{
    list.rc <- list()
    for (i in seq_len(length(genTable$filename))){
      list.rc[[i]]<- rcRun(genTable$filename[i],cytosine, posteriorMaxFilter, genTable)
    }
  }
  RCsaveResult(list.rc,cytosine,posteriorMaxFilter)
  cat(stopTime(mt))

}


rcRun <- function(filename,cytosine, posteriorMaxFilter, genTable){
  file <- RC.dataRead(filename,cytosine, posteriorMaxFilter)
  name <-  getNames(filename,genTable)
  mean.rc <-floorDec(as.numeric(mean(file$rc.meth.lvl)),5)
  res<-list(name,mean.rc)
  return(res)
}


RCsaveResult<-function(list.rc,cytosine,posteriorMaxFilter){
  tmp_dmr <- data.frame(matrix(ncol = 3, nrow =1 ))
  x <- c("Sample_name", "context", "rc.meth.lvls")
  colnames(tmp_dmr) <- x
  mainRC<-NULL
  for (nlist in seq_len(length(list.rc))){
    tmp_dmr$Sample_name<-list.rc[[nlist]][[1]]
    tmp_dmr$context<-cytosine
    tmp_dmr$rc.meth.lvls <-list.rc[[nlist]][[2]]
    mainRC<-rbind(mainRC,tmp_dmr)
  }
  mainRC<-mainRC[mixedorder(mainRC$Sample_name),]
  saveFile <- paste0(getwd(), "/", "AB-methProp-", cytosine, "-", posteriorMaxFilter, ".csv")
  fwrite(mainRC,file = saveFile ,quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
  cat(paste0("Methylation proportions results saved in: ",saveFile,"\n"))

}


