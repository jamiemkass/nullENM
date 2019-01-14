#' @export
#'

mod.args <- list()
x <- nullSDMs(occs, envs, bg, occs.folds, bg.folds, envs.folds, dismo::maxent, list(fc = "LQ", rm = 1), 10)
x <- nullSDMs(occs, envs, bg, occs.folds, bg.folds, envs.folds, maxnet::maxnet, list(fc = "LQ", rm = 1), 10)

nullSDMs <- function(occs, envs, bg, occs.folds, bg.folds, envs.folds, mod.fun, mod.args, no.iter) {

  # define new raster (envs.cv) with values equal to folds based on envs.folds
  envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
  envs.cellNo <- raster::extract(envs[[1]], envs.xy, cellnumbers = TRUE)
  envs.cv <- envs[[1]]
  envs.cv[envs.cellNo[,1]] <- envs.folds

  # get number of occurrence points by fold
  occs.folds.tbl <- table(occs.folds)

  # get total number of folds (k)
  nk <- length(unique(occs.folds))

  # initialize list to hold stats for each null model
  null.iters <- list()

  # set up empty vectors for stats
  stats.k <- data.frame(auc.test = numeric(nk), auc.diff = numeric(nk),
                        or.min = numeric(nk), or.10 = numeric(nk))

  # get environmental values for occs and bg
  occs.vals <- raster::extract(envs, occs)
  bg.vals <- raster::extract(envs, bg)

  # get model function name
  mod.fun.name <- as.character(substitute(mod.fun))[3]

  if(mod.fun.name == "maxnet") {
    mod.args$data <- as.data.frame(rbind(occs.vals, bg.vals))
    mod.args$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    mod.args$f <- maxnet::maxnet.formula(p, x, classes = tolower(mod.args$fc))
    mod.args$regmult <- mod.args$rm
    mod.args$rm <- NULL
    mod.args$fc <- NULL
    # m <- do.call(mod.fun, mod.args)
  }

  if(mod.fun.name == "maxent") {
    mod.args$args <- ENMeval::make.args(mod.args$rm, mod.args$fc)[[1]]
    mod.args$x <- as.data.frame(rbind(occs.vals, bg.vals))
    mod.args$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    mod.args$rm <- NULL
    mod.args$fc <- NULL
    # m <- do.call(mod.fun, mod.args)
  }

  # iteratively build null models
  for(i in 1:no.iter) {

    for(k in 1:nk) {
      envs.k <- envs.cv
      envs.k[envs.k == k] <- NA
      occs.null <- dismo::randomPoints(envs.k, sum(occs.folds.tbl[-k]))
      bg.train <- bg.vals[bg.folds != k, ]
      occs.test <- occs[occs.folds == k, ]

    }
  }

  return()

}

dismo::maxent(x, p, args = c(args.i, userArgs),
              factors = categoricals)
