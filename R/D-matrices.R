#' Constructing D-Matrices
#'
#' Estimating epimutation rates from high-throughput DNA methylation data
#'
#' @param genTable Generation table name, you can find sample file in
#' "extdata" called "generations.fn"
#' @param cytosine Type of cytosine (CHH/CHG/CG)
#' @param posteriorMaxFilter Filter value, based on posteriorMax
#' ex: >= 0.95 or 0.99
#' @import dplyr
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom  data.table as.data.table
#' @importFrom  stringr str_replace_all
#' @importFrom  gtools mixedorder
#' @importFrom  utils combn
#' @return generating divergence matrices file.
#' @export
#' @examples
#'## Get some toy data
#' file <- system.file("extdata","generations.fn", package="AlphaBeta")
#' df<-read.csv(file)
#' df$filename<-sub("^",paste0(dirname(file),"/"),df$filename )
#' write.csv(df, file = paste0(dirname(file),"/tm_generations.fn"),row.names=FALSE,quote=FALSE)
#' file <- system.file("extdata","tm_generations.fn", package="AlphaBeta")
#' dMatrix(file, "CG", 0.99)


dMatrix <- function(genTable, cytosine, posteriorMaxFilter) {
  # checking errors
  inputCheck(genTable, cytosine, posteriorMaxFilter)
  genTable <- fread(genTable)
  mt <- startTime("Preparing data-sets...\n")
  pairs <- combn(genTable$filename, 2)
  final_ds <- runMatrix(pairs, cytosine, posteriorMaxFilter, genTable)
  dMsaveResult(final_ds,cytosine,posteriorMaxFilter)
  cat("Generating d-matrics done.\n")
  cat(stopTime(mt))
 # rm(list=ls())
}


runMatrix <- function(pairs, cytosine, posteriorMaxFilter,genTable){
 flag=TRUE
 pair_len <- length(pairs)/2
 for (i in seq_len(pair_len)){
  df <-pairs[,i]
  name_ds<- getNames( df[1],genTable)
  name_ds[2]<- getNames(df[2],genTable)
  cat(paste0("Running: ",name_ds[[1]], " and ",
   name_ds[[2]], " ( ", i , " out of ",length(pairs)/2," pairs )"),"\n")

  cytosine<-as.character(cytosine)
  file_A <-DM.dataRead(df[1],cytosine, posteriorMaxFilter)
  file_B <-DM.dataRead(df[2],cytosine, posteriorMaxFilter)
  # check ad replace pattern in Data-set A
  returnFile <- statusStringCheck(file_A, file_B)
  file_A <- returnFile[[1]]
  file_B <- returnFile[[2]]
  cat("Computing divergence matrix...\n")
  file_A$seqnames<-as.character(file_A$seqnames)
  file_B$seqnames<-as.character(file_B$seqnames)

  tmp_db <- as.data.table(
            inner_join(file_A, file_B, by = c("seqnames","start","strand")))
  #set status 0=rows is same, 1=M/U 2=I
  rm(file_A,file_B,returnFile)

  tmp_db$state <- ifelse(tmp_db$status.x==tmp_db$status.y,0,
                  (ifelse((tmp_db$status.x=="I" | tmp_db$status.y=="I" ),2,1)))

  # substract number of Intermediate in data-set
  number_none_inter <- sum(tmp_db$state==1)     #M --> U OR U --> M
  number_intermediate <- sum(tmp_db$state==2)   #Intermediate
  Total <- NROW(tmp_db$state)                   #Total rows
  # number of (M --> U or U --> M  and 1/2 I ) / total of rows
  D <- (number_none_inter+(number_intermediate*0.5))/Total
  # reformat
  D <-floorDec(as.numeric(D),5)
    if (flag == TRUE){
       tmp_big<- data.frame(matrix(ncol = 3, nrow = 1))
       tmp_big$X1<-name_ds[[1]]
       tmp_big$X2<-name_ds[[2]]
       tmp_big$X3<-D
       flag = FALSE
    } else {
    tmp<-NULL
    tmp<-list(name_ds[[1]],name_ds[[2]],D)
    tmp_big<-rbind(tmp_big,tmp)
    rm(tmp_db)
    }
   cat(paste0(name_ds[[1]], " and ",name_ds[[2]], " is done! \n"))
   cat("|--------------------------------------------------|\n")
   name_ds<-NULL
 }
    return(tmp_big)
}


dMsaveResult<-function(final_ds,cytosine,posteriorMaxFilter){

  cat("Writing to the files:\n")
  final_ds<-final_ds[mixedorder(final_ds$X1),]
  colnames(final_ds)<-(c("pair.1", "pair.2", "D.value"))
  saved_file <- paste0(getwd(), "/", "AB-dMatrix-", cytosine, "-", posteriorMaxFilter, ".csv")
  fwrite(final_ds, file = saved_file , quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  cat(paste0("Divergence values saved in: ",saved_file, "\n"))


}

