#---------------------------------------------
# DATA Preparing - reading file
#---------------------------------------------



DM.dataRead <- function(dFile, cytosine, posteriorMaxFilter) {
    command <- sprintf(paste0("grep --text -w ", cytosine , " %s"),
      paste(dFile, collapse = ""))

    if (.Platform$OS.type == "unix") {
      tmp <-
        fread(cmd = command,
          select = c("V1", "V2", "V3", "V4", "V7", "V8"))
      cnames <-
        c("seqnames",
          "start",
          "strand",
          "context",
          "posteriorMax",
          "status")
      colnames(tmp) <- cnames

    } else{
      tmp <- fread(
        dFile,
        skip = 0,
        sep = '\t',
        select = c(
          "seqnames",
          "start",
          "strand",
          "context",
          "posteriorMax",
          "status"
        ),
        showProgress = TRUE
      )

      tmp <- tmp %>% filter(tmp$context == cytosine)
    }

    tmp <-
      tmp %>% filter(tmp$posteriorMax >= posteriorMaxFilter &
          tmp$seqnames != "M" &  tmp$seqnames != "C")


    drops <- c("context", "posteriorMax")
    tmp <- tmp %>% select(!(drops))
    #tmp <- tmp[ !(names(tmp) %in% drops)]
    # filter out based on context & posteriorMax

    return(tmp)
}


RC.dataRead <- function(dFile, cytosine, posteriorMaxFilter) {
    command <- sprintf(paste0("grep --text -w ", cytosine , " %s"),
      paste(dFile, collapse = ""))

    if (.Platform$OS.type == "unix") {
      tmp <- fread(cmd = command, select = c("V1", "V4", "V7", "V9"))
      cnames <-
        c("seqnames", "context", "posteriorMax", "rc.meth.lvl")
      colnames(tmp) <- cnames

    } else{
      tmp <- fread(
        dFile,
        skip = 0,
        sep = '\t',
        select = c("seqnames", "context", "posteriorMax", "rc.meth.lvl"),
        showProgress = TRUE
      )

      tmp <- tmp %>% filter(tmp$context == cytosine)
    }

    tmp <-
      tmp %>% filter(tmp$posteriorMax >= posteriorMaxFilter &
          tmp$seqnames != "M" &  tmp$seqnames != "C")

    # filter out based on context & posteriorMax
    drops <- c("context", "posteriorMax")
    #tmp <- tmp[,!(names(tmp) %in% drops)]
    tmp <- tmp %>% select(!(drops))
    return(tmp)
}

# take 5 digit of decimal value posteriorMax column
floorDec<-function(rc.Mean ,x){
  y <- function(x, level=1) round(x - 5*10^(-level-1), level)
  res <-y(as.numeric(rc.Mean),x)
  return(res)
}
