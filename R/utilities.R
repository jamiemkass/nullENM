# get total number of parameters
no.params <- function(mod, mod.name) {
  if(mod.name == 'maxnet') {
    return(length(mod$betas))
  }else if(mod.name == "maxent.jar") {
    lambdas <- mod@lambdas[1:(length(mod@lambdas)-4)]
    countNonZeroParams <- function(x) if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
    return(sum(unlist(sapply(lambdas, countNonZeroParams))))
  }
}

timeCheck <- function(start.time) {
  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  mins <- ifelse(t.min == 1, "minute", "minutes")
  msg <- paste(t.min, mins, round(t.sec, 1), "seconds")
  return(msg)
}

# # define a corrected variance function
# corrected.var <- function(x, nk){
#   sum((x - mean(x))^2) * ((nk-1)/nk)
# }
