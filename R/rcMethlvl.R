#' Calculating rc.Meth.lvl
#'
#' Estimating epimutation rates from high-throughput DNA methylation data
#' @param genTable Generation table name, you can find sample file in
#' "extdata" called "generations.fn"
#' @param cytosine Type of cytosine (CHH/CHG/CG)
#' @param posteriorMaxFilter Filter value, based on posteriorMax
#' @import      dplyr
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom  data.table as.data.table
#' @importFrom  stringr str_replace_all
#' @importFrom  gtools mixedorder
#' @importFrom  BiocParallel bplapply
#' @importFrom  BiocParallel SnowParam
#' @return rc meth lvl.
#' @export
#' @examples
#'## Get some toy data
#' file <- system.file("extdata","tm_generations.fn", package="AlphaBeta")
#' rc.meth.lvl(file, "CG", 0.99)


# running join
rc.meth.lvl <- function(genTable, cytosine, posteriorMaxFilter){

      inputCheck(genTable, cytosine, posteriorMaxFilter)
      genTable <- fread(genTable)
      mt <- startTime("Preparing data-sets...\n")

      list.rc<-bplapply(genTable$filename,cytosine=cytosine,posteriorMaxFilter= posteriorMaxFilter,
        genTable= genTable, rcRun, BPPARAM = SnowParam(exportglobals = FALSE))

      RCsaveResult(list.rc,cytosine,posteriorMaxFilter)
      cat(stopTime(mt))

}


rcRun <- function(filename, cytosine, posteriorMaxFilter, genTable){
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


