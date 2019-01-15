#' @export
#'

mod.args <- list(fc = "LQ", rm = 1)
# x <- nullSDMs(occs, envs, bg, occs.folds, bg.folds, envs.folds, dismo::maxent, mod.args, 10)
# x <- nullSDMs(occs, envs, bg, occs.folds, bg.folds, envs.folds, maxnet::maxnet, list(fc = "LQ", rm = 1), 10)

nullSDMs <- function(occs, envs, bg, occs.folds, bg.folds, envs.folds, mod.fun, mod.args, no.iter, categoricals) {

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

  # get environmental values for occs and bg
  occs.vals <- as.data.frame(raster::extract(envs, occs))
  bg.vals <- as.data.frame(raster::extract(envs, bg))

  # get model function name
  mod.fun.name <- as.character(substitute(mod.fun))[3]



  # iteratively build null models
  for(i in 1:no.iter) {
    # set up empty vectors for stats
    stats.i <- data.frame(auc.test = numeric(nk), auc.diff = numeric(nk),
                          or.min = numeric(nk), or.10 = numeric(nk))

    # sample random occurrences
    occs.null <- list()
    for(k in 1:nk) {
      envs.k <- envs.cv
      envs.k[envs.k != k] <- NA
      occs.null.xy <- dismo::randomPoints(envs.k, occs.folds.tbl[k])
      occs.null[[k]] <- as.data.frame(extract(envs, occs.null.xy))
    }
    occs.null.all <- do.call(rbind, occs.null)


    mod.args.i <- buildModArgs(mod.fun.name, mod.args, occs.null.all, bg.vals)
    mod.i <- do.call(mod.fun, mod.args.i)
    auc.train.i <- dismo::evaluate(occs.null.all, bg.vals, mod.i)@auc

    for(k in 1:nk) {

      occs.null.train <- occs.null[[k]]
      bg.train <- bg.vals[bg.folds != k, ]
      occs.test <- occs.vals[occs.folds == k, ]

      # convert fields for categorical data to factor class
      if(!is.null(categoricals)) {
        for (i in 1:length(categoricals)) {
          occs.null.train[, categoricals[i]] <- as.factor(occs.null.train[, categoricals[i]])
          bg.train[, categoricals[i]] <- as.factor(bg.train[, categoricals[i]])
        }
      }

      # build custom mod.args for this iteration
      mod.args.k <- buildModArgs(mod.fun.name, mod.args, occs.null.train, bg.train)
      # build the model
      mod.k <- do.call(mod.fun, mod.args.k)
      # evaluate on withheld data (occs.test)
      auc.test.k <- dismo::evaluate(occs.test, bg.train, mod.k)@auc

    }
  }

  return()

}

#dismo::maxent(x, p, args = c(args.i, userArgs),factors = categoricals)
