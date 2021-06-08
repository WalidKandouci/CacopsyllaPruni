
## ###########################################
## Take arguments from script, if available ##
CA <- commandArgs(TRUE)
if (length(CA)==0) {
  UseScript <- FALSE
} else {
  UseScript <- TRUE
}

if (UseScript) {
  print(CA)
  print(iModel <- as.integer(CA)[1])
  print(qsubID <- as.integer(CA)[2])
} else {
  qsubID <- 123456789
}
