inputCheck <- function(...) {
  # getting var parameters
  var <- list(...)
  # input file checking

  # check if genTable is exist then check file exist one-by-one
  if (!file.exists(var[[1]])){
    stop("Sample file dosen't exist.")
  }
  gen_tbl <- fread(var[[1]])
  gen_tbl <- gen_tbl %>% filter(gen_tbl$meth=="Y")
  for (i in seq_len(NROW(gen_tbl))){
    if (!file.exists(gen_tbl[[i,1]])){
      stop(paste0("File ",gen_tbl[[i,1]], " Not Exist!"))
    }
  }

  # cytosines part checking
  list_cyt <- c("CHH", "CHG", "CG")
  cytosines <- var[[2]]
  if (!cytosines %in% list_cyt | length(cytosines) > 1) {
    stop("Please enter the valid Cytosines. (CHH/CHG/CG)")
  }

  # posteriorMaxFiltervalue part checking
  pMFilter <- as.numeric(var[[3]])
  if (pMFilter <= 0.80 | pMFilter > 1) {
    stop("posteriorMax value must be < 1 and >= 0.80 ")
  }
  return(0)
}

statusStringCheck <-  function(file_A, file_B){
  list_status <- c("Unmethylated", "Intermediate", "Methylated")

  strTocheckFileA <- utils::head(file_A$status[1])
  if (strTocheckFileA %in% list_status) {
    file_A$status <- str_replace_all(file_A$status,
                      pattern = "Unmethylated", replacement = "U")
    file_A$status <- str_replace_all(file_A$status,
                      pattern = "Intermediate", replacement = "I")
    file_A$status <- str_replace_all(file_A$status,
                      pattern = "Methylated", replacement = "M")
  }

  strTocheckFileB <- utils::head(file_B$status[1])
  if (strTocheckFileB %in% list_status) {
    # replace pattern in Data-set B
    file_B$status <- str_replace_all(file_B$status,
                      pattern = "Unmethylated", replacement = "U")
    file_B$status <- str_replace_all(file_B$status,
                      pattern = "Intermediate", replacement = "I")
    file_B$status <- str_replace_all(file_B$status,
                      pattern = "Methylated", replacement = "M")
  }
  return(list(file_A, file_B))

}

getNames <- function(nameDF,genTable){
  #genTable <- fread(genTable)
  strSearch<-nameDF
  tmp_A <- genTable[genTable$filename == strSearch,][,2]
  #tmp_B <- genTable[genTable$filename == strSearch,][,3]
  #if (tmp_B == "" | is.na(tmp_B) | is.null(tmp_B) ){
  #  tmp_B<-""
  #  gen_name <- paste0(tmp_A,tmp_B)
  #}else{
  #  gen_name <- paste0(tmp_A,"-",tmp_B)
  #}
  #return(gen_name)
  return(tmp_A)
}
