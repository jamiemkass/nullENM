#' @export

buildModArgs <- function(mod.fun.name, mod.args, occs.vals, bg.vals) {
  if(mod.fun.name == "maxnet") {
    mod.args$data <- rbind(occs.vals, bg.vals)
    mod.args$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    mod.args$f <- maxnet::maxnet.formula(mod.args$p, mod.args$data, classes = tolower(mod.args$fc))
    mod.args$regmult <- mod.args$rm
    mod.args$rm <- NULL
    mod.args$fc <- NULL
  }

  if(mod.fun.name == "maxent") {
    mod.args$args <- ENMeval::make.args(mod.args$rm, mod.args$fc)[[1]]
    mod.args$x <- as.data.frame(rbind(occs.vals, bg.vals))
    mod.args$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    mod.args$rm <- NULL
    mod.args$fc <- NULL
  }

  return(mod.args)
}
