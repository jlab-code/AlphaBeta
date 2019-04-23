#' Constructing D-Matrices
#'
#' Estimating epimutation rates from high-throughput DNA methylation data
#'
#' @param genTable Generation table name, you can find sample file in
#' "extdata" called "generations.fn"
#' @param cytosines Type of cytosines (CHH/CHG/CG)
#' @param posteriorMaxFilter Filter value, based on posteriorMax
#' ex: >= 0.95 or 0.99
#' @import dplyr
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom  data.table as.data.table
#' @importFrom  stringr str_replace_all
#' @return generating divergence matrices file.
#' @export
#' @examples
#'## Get some toy data
#' file <- system.file("extdata","generations.fn", package="alphabeta")
#'
#' dMatrices(file, "CG", 0.99)


dMatrices <- function(genTable, cytosines, posteriorMaxFilter) {
  # checking errors

  inputCheck(genTable, cytosines, posteriorMaxFilter)

  gen_tbl <- fread(genTable)

  mt <- startTime("Preparing data-sets...\n")
  pairs <- utils::combn(gen_tbl$samplename, 2)

  final_ds <- runMatrix(pairs, cytosines, posteriorMaxFilter, genTable)
  cat("Generating d-matrics done.\n")
  cat("Writing to the files:\n")
  colnames(final_ds)<-(c("pair-1", "pair-2", "D-value", "cytosines", "filter"))
  saved_file <-
    paste0(getwd(), "/", "d-matrices-", cytosines, "-", posteriorMaxFilter, ".csv")
  fwrite(
    final_ds,
    file = saved_file ,
    quote = FALSE,
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE
  )
  cat(paste0(saved_file, "\n"))
  cat(stopTime(mt))
}

# -------------------------------------------------------------
# getting file names
#-------------------------------------------------------------
getNames <- function(genTable, nameDF){
    gen_tbl <- fread(genTable)
    strSearch<-nameDF
    tmp_A <- gen_tbl[gen_tbl$samplename == strSearch,][,2]
    tmp_B <- gen_tbl[gen_tbl$samplename == strSearch,][,3]
    if (tmp_B == "" | is.na(tmp_B) | is.null(tmp_B) ){
    tmp_B<-""
    gen_name <- paste0(tmp_A,tmp_B)
    }else{
    gen_name <- paste0(tmp_A,"-",tmp_B)
    }
    return(gen_name)
}

runMatrix <- function(pairs, cytosines, posteriorMaxFilter, genTable){
 flag=TRUE
 pair_len <- length(pairs)/2
 for (i in seq_len(pair_len)){
  df <-pairs[,i]
  name_ds<- getNames(genTable, df[1])
  name_ds[2]<- getNames(genTable, df[2])
  cat(paste0("Running: ",name_ds[[1]], " and ",
   name_ds[[2]], " ( ", i , " out of ",length(pairs)/2," pairs )"),"\n")

  cat("Reading data-set and constructing data-frames ...\n")
  file_A <- fread(df[1], skip = 0, sep = '\t',
       select=c("seqnames","start","strand","context","posteriorMax","status"),
       showProgress=FALSE)
  file_B <- fread(df[2], skip = 0, sep = '\t',
       select=c("seqnames","start","strand","context","posteriorMax","status"),
       showProgress=FALSE)

  # filter out based on context & posteriorMax
  file_A <- file_A %>%
  filter(file_A$context == cytosines & file_A$posteriorMax >= posteriorMaxFilter )
  file_B <- file_B %>%
  filter(file_B$context == cytosines & file_B$posteriorMax >= posteriorMaxFilter )
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
  floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
  D <-floor_dec(as.numeric(D),6)
    if (flag == TRUE){
       tmp_big<- data.frame(matrix(ncol = 5, nrow = 1))
       tmp_big$X1<-name_ds[[1]]
       tmp_big$X2<-name_ds[[2]]
       tmp_big$X3<-D
       tmp_big$X4<-cytosines
       tmp_big$X5<-posteriorMaxFilter
       flag = FALSE
    } else {
    tmp<-NULL
    tmp<-list(name_ds[[1]],name_ds[[2]],D,cytosines,posteriorMaxFilter)
    tmp_big<-rbind(tmp_big,tmp)
    rm(tmp_db)
    }
   cat(paste0(name_ds[[1]], " and ",name_ds[[2]], " is done! \n"))
   cat("|--------------------------------------------------|\n")
   name_ds<-NULL
 }
    return(tmp_big)
}




