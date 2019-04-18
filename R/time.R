startTime <- function(...) {
  x <- paste0(..., collapse='')
  message(x, appendLF=FALSE)
  st <- proc.time()
  return(st)
}
stopTime <- function(st) {
  sttime <- proc.time() - st
  message("Total time: ", round(sttime[3],2), "s")
}
