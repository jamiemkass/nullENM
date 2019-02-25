#' @export
#'

# for split evaluation, label training occs "1" and independent evaluation occs "2" in partitions
nullSDMs <- function(occs, envs, bg, partitions, mod.name, mod.args, no.iter,
                     evalType = c("split", "kfold", "kspatial"),
                     categoricals = NULL, abs.auc.diff = FALSE, removeMxTemp = TRUE) {

  # record start time
  start.time <- proc.time()
  message("Beginning null SDM analysis...")

  # changing stack to brick to speed up analysis
  envs <- raster::brick(envs)

  # if split evaluation, partition number equals 1
  if(evalType == "split") {
    nk <- 1
  }else{
    # get total number of partitions (k)
    nk <- length(unique(partitions$occ.grp))
  }
  # if spatial k-fold evaluation, spatially partition envs
  if(evalType == "kspatial") {
    # define new raster (envs.cv) with values equal to partitions based on partitions$envs.grp
    envs.xy <- raster::rasterToPoints(envs[[1]], spatial = TRUE)
    envs.cellNo <- raster::extract(envs[[1]], envs.xy, cellnumbers = TRUE)
    envs.cv <- envs[[1]]
    envs.cv[envs.cellNo[,1]] <- partitions$envs.grp
  }

  # get number of occurrence points by partition
  occs.parts.tbl <- table(partitions$occ.grp)

  # get environmental values for occs and bg
  occs.vals <- as.data.frame(raster::extract(envs, occs))
  bg.vals <- as.data.frame(raster::extract(envs, bg))

  # convert fields for categorical data to factor class
  if(!is.null(categoricals)) {
    for(j in 1:length(categoricals)) {
      occs.vals[, categoricals[j]] <- as.factor(occs.vals[, categoricals[j]])
      bg.vals[, categoricals[j]] <- as.factor(bg.vals[, categoricals[j]])
    }
  }

  # get model function name (only Maxent functions available now)
  if(mod.name == "maxent.jar") {
    mod.fun <- dismo::maxent
    # create temp directory to store maxent.jar output, for potential removal later
    tmpdir <- paste(tempdir(), runif(1,0,1), sep = "/")
    dir.create(tmpdir, showWarnings = TRUE, recursive = FALSE)
  }else if(mod.name == "maxnet") {
    mod.fun <- maxnet:maxnet

  }else{
    message('Only Maxent functions available now. Please choose either "maxent.jar" or "maxnet".')
    return()
  }

  # initialize data frames to collect evaluation stats for each partition
  # per iteration and their averages
  all.cnames <- c("auc.train", "auc.test", "auc.diff", "or.min", "or.10")
  k.cnames <- c("auc.test", "auc.diff", "or.min", "or.10")
  null.cnames <- c("auc.train", "mean.auc.test", "sd.auc.test",
                         "mean.auc.diff", "sd.auc.diff", "mean.or.min", "sd.or.min",
                         "mean.or.10", "sd.or.10", "nparam")
  all.rnames <- c("real.mean", "real.sd", "null.mean", "null.sd", "zscore", "pvalue")
  all.stats <- data.frame(matrix(nrow = 6, ncol = length(all.cnames),
                                 dimnames = list(all.rnames, all.cnames)))
  # real.stats <- data.frame(matrix(nrow = 1, ncol = 2,
                                  # dimnames = list(NULL, c("auc.train", "nparam"))))
  null.stats <- data.frame(matrix(nrow = no.iter, ncol = length(null.cnames),
                                  dimnames = list(NULL, null.cnames)))

  ############################## #
  # build real model ####
  ############################## #

  mod.args.real <- buildModArgs(mod.name, mod.args, occs.vals, bg.vals)
  mod.real <- do.call(mod.fun, mod.args.real)
  # calculate training auc
  all.stats["real.mean", "auc.train"] <- dismo::evaluate(occs.vals, bg.vals, mod.real)@auc
  # real.stats$nparam <- no.params(mod.real, mod.name)
  kstats <- data.frame(matrix(nrow = nk, ncol = length(k.cnames),
                              dimnames = list(NULL, k.cnames)))

  for(k in 1:nk) {
    if(evalType == "kfold" | evalType == "kspatial") {
      occs.train.real.k <- occs.vals[partitions$occ.grp != k, ]
      occs.test.real.k <- occs.vals[partitions$occ.grp == k, ]
      bg.train.k <- bg.vals[partitions$bg.grp != k, ]
    }else if(evalType == "split") {
      occs.train.real.k <- occs.vals[partitions$occ.grp != 2,]
      occs.test.real.k <- occs.vals[partitions$occ.grp != 1,]
      bg.train.k <- bg.vals
    }
    mod.args.real.k <- buildModArgs(mod.name, mod.args, occs.train.real.k, bg.train.k)
    mod.real.k <- do.call(mod.fun, mod.args.real.k)
    kstats[k,] <- evalStats(occs.train.real.k, bg.train.k, occs.test.real.k,
                               mod.real.k, abs.auc.diff)
    message(sprintf("Completed real model and evaluation for partition %i.", k))
  }
  # fill in rest of real model statistics
  all.stats[1, 2:5] <- apply(kstats, 2, mean)
  all.stats[2, 2:5] <- apply(kstats, 2, sd)

  ############################## #
  # build null models ####
  ############################## #

  # initialize list to record stats for null iteration i
  null.stats.iters <- list()

  for(i in 1:no.iter) {
    # initialize data frame for partition statistics for current iteration
    null.stats.iters[[i]] <- data.frame(matrix(nrow = nk, ncol = length(k.cnames),
                                           dimnames = list(NULL, k.cnames)))

    # sample random occurrences
    occs.null <- list()
    if(evalType == "kfold") {
      # if kfold evaluation, randomly sample the same number of training
      # occs over each k partition of envs
      for(k in 1:nk) {
        occs.null.xy <- dismo::randomPoints(envs, occs.parts.tbl[k])
        occs.null[[k]] <- as.data.frame(raster::extract(envs, occs.null.xy))
      }
    }else if(evalType == "kspatial") {
      # if kfold evaluation, randomly sample the same number of training
      # occs over each k partition of envs
      for(k in 1:nk) {
        envs.k <- envs.cv
        envs.k[envs.k != k] <- NA
        occs.null.xy <- dismo::randomPoints(envs.k, occs.parts.tbl[k])
        occs.null[[k]] <- as.data.frame(raster::extract(envs, occs.null.xy))
      }
    }else if(evalType == "split") {
      # if split evaluation, randomly sample the same number of training
      # occs over the envs extent only once
      occs.null.xy <- dismo::randomPoints(envs, occs.parts.tbl[1])
      occs.null[[1]] <- as.data.frame(raster::extract(envs, occs.null.xy))
    }

    # convert fields for categorical data to factor class
    if(!is.null(categoricals)) {
      for(k in 1:nk) {
        for(j in 1:length(categoricals)) {
          occs.null[[k]][, categoricals[j]] <- as.factor(occs.null[[k]][, categoricals[j]])
        }
      }
    }

    occs.null.all <- do.call(rbind, occs.null)

    mod.args.i <- buildModArgs(mod.name, mod.args, occs.null.all, bg.vals)
    mod.i <- do.call(mod.fun, mod.args.i)
    null.stats[i,]$auc.train <- dismo::evaluate(occs.null.all, bg.vals, mod.i)@auc
    null.stats[i,]$nparam <- no.params(mod.i, mod.name)

    # get training and testing data for current iteration
    for(k in 1:nk) {
      if(evalType == "kfold" | evalType == "kspatial") {
        # if kfold evaluation, use non-k training and k testing partitions
        occs.null.train <- do.call(rbind, occs.null[-k])
        bg.train <- bg.vals[partitions$bg.grp != k, ]
        occs.test <- occs.vals[partitions$occ.grp == k, ]
      }else if(evalType == "split") {
        # if split evaluation, use single random occs dataset and real
        # independent testing dataset
        occs.null.train <- occs.null[[1]]
        bg.train <- bg.vals
        occs.test <- occs.vals[partitions$occ.grp != 1,]
      }

      # build custom mod.args for this iteration
      mod.args.k <- buildModArgs(mod.name, mod.args, occs.null.train, bg.train)
      # build the model
      mod.k <- do.call(mod.fun, mod.args.k)

      null.stats.iters[[i]][k,] <- evalStats(occs.null.train, bg.train, occs.test,
                                         mod.k, abs.auc.diff)

      message(sprintf("Completed null model and evaluation for partition %i of iteration %i.", k, i))
    }
    # average partition statistics
    # for split partition, sd will be NA because input is single fold statistic
    null.stats.means <- apply(null.stats.iters[[i]], 2, mean)
    null.stats.sds <- apply(null.stats.iters[[i]], 2, sd)
    # alternate the above values for data frame and assign to table
    null.stats[i,2:9] <- c(rbind(null.stats.means, null.stats.sds))
  }

  # calculate means and standard deviations of mean k-fold values for null stats
  all.stats["null.mean",] <- apply(null.stats[c(1,2,4,6,8)], 2, mean)
  all.stats["null.sd",] <- apply(null.stats[c(1,2,4,6,8)], 2, sd)
  all.stats["zscore",] <- (all.stats["real.mean",] - all.stats["null.mean",]) / all.stats["null.sd",]
  all.stats["pvalue", 1:2] <- 1 - sapply(all.stats["zscore", 1:2], pnorm)
  all.stats["pvalue", 3:5] <- sapply(all.stats["zscore", 3:5], pnorm)

  out <- list(all.stats = all.stats, null.stats = null.stats, null.stats.iters = null.stats.iters)

  # optionally remove temp directory for maxent.jar
  if(mod.name == "maxent.jar" & removeMxTemp == TRUE) unlink(tmpdir, recursive = T, force = T)

  timed <- proc.time() - start.time
  t.min <- floor(timed[3] / 60)
  t.sec <- timed[3] - (t.min * 60)
  message(paste("Null SDM analysis completed in", t.min, "minutes", round(t.sec, 1), "seconds."))

  return(out)
}
