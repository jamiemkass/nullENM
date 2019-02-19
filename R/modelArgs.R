#' @export

#' @export
model.args <- function(tune.tbl.i, mod.name, occs.vals, bg.vals, other.args) {

  out <- list()
  # define data
  d <- rbind(occs.vals, bg.vals)
  # define response
  p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  # maxent.jar
  if(mod.name == "maxent") {
    out$x <- d
    out$p <- p
    out$args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
    if(!grepl("L", tune.tbl.i$fc)) out$args <- c(out$args, "nolinear")
    if(!grepl("Q", tune.tbl.i$fc)) out$args <- c(out$args, "noquadratic")
    if(!grepl("H", tune.tbl.i$fc)) out$args <- c(out$args, "nohinge")
    if(!grepl("P", tune.tbl.i$fc)) out$args <- c(out$args, "noproduct")
    if(!grepl("T", tune.tbl.i$fc)) out$args <- c(out$args, "nothreshold")
    out$args <- c(out$args, paste0("betamultiplier=", tune.tbl.i$rm, sep=""))
  }
  # maxnet
  if(mod.name == "maxnet") {
    out$data <- d
    out$p <- p
    out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.tbl.i$fc))
    out$regmult <- tune.tbl.i$rm
  }

  # add other args
  out <- c(out, other.args)

  return(out)
}
