#' @export
model.args <- function(mod.name, mod.args, occs.vals, bg.vals, other.args) {
  out <- list()
  # define data
  d <- rbind(occs.vals, bg.vals)
  # define response
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  # maxent.jar
  if(mod.name == "maxent.jar") {
    out$x <- d
    out$p <- p
    out$args <- c("noaddsamplestobackground", "noremoveDuplicates")
    if(!is.null(mod.args)) {
      out$args <- c(out$args, "noautofeature")
      if(!grepl("L", mod.args$fc)) out$args <- c(out$args, "nolinear")
      if(!grepl("Q", mod.args$fc)) out$args <- c(out$args, "noquadratic")
      if(!grepl("H", mod.args$fc)) out$args <- c(out$args, "nohinge")
      if(!grepl("P", mod.args$fc)) out$args <- c(out$args, "noproduct")
      if(!grepl("T", mod.args$fc)) out$args <- c(out$args, "nothreshold")
      out$args <- c(out$args, paste0("betamultiplier=", mod.args$rm, sep=""))
    }
  }
  # maxnet
  if(mod.name == "maxnet") {
    out$data <- d
    out$p <- p
    if(!is.null(mod.args)) {
      out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(mod.args$fc))
      out$regmult <- mod.args$rm
    }else{
      out$f <- maxnet::maxnet.formula(out$p, out$data)
    }
  }

  # add other args
  out <- c(out, other.args)

  return(out)
}
