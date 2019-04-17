#' Constructing D-Matrices
#'
#' Estimating epimutation rates from high-throughput DNA methylation data
#'
#' @param directory Path to the data-source directory
#' @param cytosines Type of cytosines (CHH/CHG/CG)
#' @param posteriorMaxFilter Filter value, based on posteriorMax ex: >= 0.95 or 0.99
#' @param genTable Generation table name, you can find sample file in "extdata" generations.fn
#' @import dplyr
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom data.table as.data.table
#' @importFrom  stringr str_replace_all
#' @return generating divergence matrices file.
#' @export
#' @examples
#' dMatrices("data-source directory", "CHH", 0.99, "table.fn")
#' generation table name should be in this format:
#' system.file("extdata","generations.fn", package = "alphabeta")


dMatrices<-function(directory, cytosines, posteriorMaxFilter, genTable ){
  # checking errors
  if (!dir.exists(directory)){
    stop("Directory of data-source not exists!")
  }
  # check if cytosines is correct

  list_cyt<-c("CHH","CHG","CG")
  if (!cytosines %in% list_cyt | length(cytosines) > 1 ){
    stop("Please enter the valid Cytosines. (CHH/CHG/CG)")
  }
  # check if posteriorMaxFiltervalue is correct
  pMFilter <- as.numeric(posteriorMaxFilter)
  if (pMFilter <= 0 | pMFilter > 1 ){
    stop("posteriorMax value must be < 1 and >= 0. ")
  }
  data.names <- list.files(paste0(dirname(directory),"/",basename(directory)),pattern = "*.txt|*.tsv", full.names = TRUE)
  gen_tbl <- fread(genTable)
  if (length(data.names)!= NROW(gen_tbl)){
    stop("Number of file is not matching with data file generation names. ")
  }
  mt <- startTime("Preparing data-sets...\n")
  pairs<-combn(data.names, 2)
  final_ds <- runMatrix(pairs, cytosines, pMFilter, genTable)
  cat("Generating d-matrics done.\n")
  cat("Writing to the files:\n")
  colnames(final_ds)<-(c("pair-1","pair-2","D-value","cytosines","filter"))
  saved_file <- paste0(getwd(),"/","d-matrices-",cytosines,"-",pMFilter,".csv")
  fwrite(final_ds,file = saved_file , quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  cat(paste0(saved_file,"\n"))
  cat(paste0("Total time: ",stopTime(mt)))
}

# -------------------------------------------------------------
# getting file names
#-------------------------------------------------------------
getNames <- function(genTable, nameDF){
  gen_tbl <- fread(genTable)
  strSearch<-basename(nameDF)
  tmp_A <- gen_tbl[gen_tbl$samplename == strSearch,][,2]
  tmp_B <- gen_tbl[gen_tbl$samplename == strSearch,][,3]
  if (tmp_B == ""){
  gen_name <- paste0(tmp_A,tmp_B)
  }else{
    gen_name <- paste0(tmp_A,"-",tmp_B)
  }
  return(gen_name)
}
runMatrix <- function(pairs, cytosines, pMFilter, genTable){

  flag=TRUE
  for (i in 1:(length(pairs)/2)){
      df <-pairs[,i]
      name_ds<- getNames(genTable, df[1])
      name_ds[2]<- getNames(genTable, df[2])
      cat(paste0("Running: ",name_ds[[1]], " and ",name_ds[[2]], " ( ", i , " out of ",length(pairs)/2," pairs )"),"\n")

      cat("Reading data-set and constructing data-frames ...\n")
      file_A <- fread(df[1], skip = 0, sep = '\t',select=c("seqnames","start","strand","context","posteriorMax","status"), showProgress=FALSE)
      file_B <- fread(df[2], skip = 0, sep = '\t',select=c("seqnames","start","strand","context","posteriorMax","status"), showProgress=FALSE)
      # filter out based on context & posteriorMax
      file_A <- file_A %>% filter(context == cytosines & posteriorMax >= pMFilter )
      file_B <- file_B %>% filter(context == cytosines & posteriorMax >= pMFilter )
      # replace pattern in Data-set A
      file_A$status <- str_replace_all(file_A$status, pattern = "Unmethylated", replacement = "U")
      file_A$status <- str_replace_all(file_A$status, pattern = "Intermediate", replacement = "I")
      file_A$status <- str_replace_all(file_A$status, pattern = "Methylated", replacement = "M")
      # replace pattern in Data-set B
      file_B$status <- str_replace_all(file_B$status, pattern = "Unmethylated", replacement = "U")
      file_B$status <- str_replace_all(file_B$status, pattern = "Intermediate", replacement = "I")
      file_B$status <- str_replace_all(file_B$status, pattern = "Methylated", replacement = "M")

      cat("Computing divergence matrix...\n")
      file_A$seqnames<-as.character(file_A$seqnames)
      file_B$seqnames<-as.character(file_B$seqnames)
      tmp_db <- as.data.table(inner_join(file_A , file_B, by = c("seqnames","start","strand")))
      ## set status 0=rows is same , 1=is Methylated/Unmethylated 2=there is Intermediate
      rm(file_A,file_B)
      tmp_db$state <- ifelse(tmp_db$status.x==tmp_db$status.y,0,(ifelse((tmp_db$status.x=="I" | tmp_db$status.y=="I" ),2,1)))

    # substract number of Intermediate in data-set
    number_none_inter<-sum(tmp_db$state==1)   # M --> U OR U --> M
    number_intermediate<-sum(tmp_db$state==2) # Intermediate
    Total=NROW(tmp_db$state)                  # Total rows
    # number of (Methylated --> Unmethylated ) / total of rows
    D=(number_none_inter+(number_intermediate*0.5))/Total
    # reformat
    floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
    D <-floor_dec(as.numeric(D),6)
    if (flag == TRUE){
      tmp_big<- data.frame(matrix(ncol = 5, nrow = 1))
      tmp_big$X1<-name_ds[[1]]
      tmp_big$X2<-name_ds[[2]]
      tmp_big$X3<-D
      tmp_big$X4<-cytosines
      tmp_big$X5<-pMFilter
      flag = FALSE
    } else {
      tmp<-NULL
      tmp<-list(name_ds[[1]],name_ds[[2]],D,cytosines,pMFilter)
      tmp_big<-rbind(tmp_big,tmp)
      rm(tmp_db)
    }
    cat(paste0(name_ds[[1]], " and ",name_ds[[2]], " is done! \n"))
    cat("|--------------------------------------------------|\n")

    name_ds<-NULL
  }
  return(tmp_big)
}




