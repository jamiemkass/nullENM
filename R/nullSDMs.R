#' @export
#'

# bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=300, has_coords=TRUE)
# occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
# occs <- occs[!duplicated(occs),]
# envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
envs <- raster::mask(envs, envs$biome)
# which(rowSums(is.na(raster::extract(envs, occs))) > 0)
# bg <- dismo::randomPoints(envs[[1]], 10000)
# mod.args <- list(fc = "LQ", rm = 1)
# folds <- ENMeval::get.block(occs, bg)
# occs.folds <- folds$occ.grp
# bg.folds <- folds$bg.grp
# envs.xy <- raster::rasterToPoints(envs[[1]], spatial = TRUE)
# envs.folds <- ENMeval::get.block(occs, envs.xy@coords)$bg.grp
# x1 <- nullSDMs(occs, envs, bg, occs.folds, bg.folds, envs.folds, dismo::maxent, mod.args, 10, "biome")
# x2 <- nullSDMs(occs, envs, bg, occs.folds, bg.folds, envs.folds, maxnet::maxnet, mod.args, 10, "biome")

nullSDMs <- function(occs, envs, bg, occs.folds, bg.folds, envs.folds,
                     mod.fun, mod.args, no.iter, categoricals=NULL, abs.auc.diff = FALSE) {

  # define new raster (envs.cv) with values equal to folds based on envs.folds
  envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
  envs.cellNo <- raster::extract(envs[[1]], envs.xy, cellnumbers = TRUE)
  envs.cv <- envs[[1]]
  envs.cv[envs.cellNo[,1]] <- envs.folds

  # get number of occurrence points by fold
  occs.folds.tbl <- table(occs.folds)

  # get total number of folds (k)
  nk <- length(unique(occs.folds))

  # get environmental values for occs and bg
  occs.vals <- as.data.frame(raster::extract(envs, occs))
  bg.vals <- as.data.frame(raster::extract(envs, bg))

  # get model function name
  mod.fun.name <- as.character(substitute(mod.fun))[3]

  # initialize data frames to collect evaluation stats for each fold
  # per iteration and their averages
  fields.avg <- c("iter", "full.auc.train", "avg.auc.train", "std.auc.train",
                  "avg.auc.test", "std.auc.test", "avg.auc.diff", "std.auc.diff",
              "avg.or.min", "std.or.min", "avg.or.10", "std.or.10")
  null.stats.avg <- data.frame(matrix(nrow = no.iter, ncol = length(fields.avg)))
  names(null.stats.avg) = fields.avg
  null.stats.i <- list()



  # iteratively build null models
  for(i in 1:no.iter) {
    # initialize data frame for fold statistics for current iteration
    fields.i <- c("auc.train", "auc.test", "auc.diff", "or.min", "or.10")
    null.stats.i[[i]] <- data.frame(matrix(nrow = nk, ncol = length(fields.i)))
    names(null.stats.i[[i]]) <- fields.i

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
      # get training and testing data for current iteration
      # (non-k training and k testing folds)
      occs.null.train <- do.call(rbind, occs.null[-k])
      bg.train <- bg.vals[bg.folds != k, ]
      occs.test <- occs.vals[occs.folds == k, ]

      # convert fields for categorical data to factor class
      if(!is.null(categoricals)) {
        for(j in 1:length(categoricals)) {
          occs.null.train[, categoricals[j]] <- as.factor(occs.null.train[, categoricals[j]])
          bg.train[, categoricals[j]] <- as.factor(bg.train[, categoricals[j]])
        }
      }

      # build custom mod.args for this iteration
      mod.args.k <- buildModArgs(mod.fun.name, mod.args, occs.null.train, bg.train)
      # build the model
      mod.k <- do.call(mod.fun, mod.args.k)

      # calculate auc on training and testing data
      auc.train.k <- dismo::evaluate(occs.null.train, bg.train, mod.k)@auc
      auc.test.k <- dismo::evaluate(occs.test, bg.train, mod.k)@auc
      # calculate auc diff
      auc.diff.k <- auc.train.k - auc.test.k
      if(abs.auc.diff == TRUE) auc.diff.k <- abs(auc.diff.k)
      # get model predictions for training and testing data
      pred.train.k <- predict(mod.k, occs.null.train)
      pred.test.k <- predict(mod.k, occs.test)
      # get 10 percentile predicted value
      occs.null.train.n <- nrow(occs.null.train)
      if(occs.null.train.n < 10) {
        pct10.train.k <- ceiling(occs.null.train.n * 0.1)
      } else {
        pct10.train.k <- floor(occs.null.train.n * 0.1)
      }
      pct10.train.k.thr <- sort(pred.train.k)[pct10.train.k]
      or10.test.k <- mean(pred.test.k < pct10.train.k.thr)
      min.train.k.thr <- min(pred.train.k)
      orMin.test.k <- mean(pred.test.k < min.train.k.thr)

      null.stats.i[[i]][k,] <- c(auc.train.k, auc.test.k, auc.diff.k,
                                 orMin.test.k, or10.test.k)

      message(sprintf("Completed model and evaluation for fold %i of iteration %i.", k, i))
    }
    # average fold statistics
    means <- apply(null.stats.i[[i]], 2, mean)
    stds <- apply(null.stats.i[[i]], 2, sd)
    # alternate the above values for data frame
    stats <- c(rbind(means, stds))
    null.stats.avg[i,] <- c(i, auc.train.i, stats)
  }

  out <- list(null.stats.avg = null.stats.avg, null.stats.i = null.stats.i)
  return(out)

}

#dismo::maxent(x, p, args = c(args.i, userArgs),factors = categoricals)
