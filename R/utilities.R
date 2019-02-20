# get total number of parameters
no.params <- function(mod, mod.name) {
  if(mod.name == 'maxnet') {
    return(length(mod$betas))
  }else if(mod.name == "maxent") {
    lambdas <- mod@lambdas[1:(length(mod@lambdas)-4)]
    countNonZeroParams <- function(x) if(strsplit(x, split=", ")[[1]][2] != '0.0') 1
    return(sum(unlist(sapply(lambdas, countNonZeroParams))))
  }
}

# # define a corrected variance function
# corrected.var <- function(x, nk){
#   sum((x - mean(x))^2) * ((nk-1)/nk)
# }
