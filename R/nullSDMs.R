#' @export
#'

# bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=300, has_coords=TRUE)
# occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
# occs <- occs[!duplicated(occs),]
# envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
# envs <- raster::mask(envs, envs$biome)
# which(rowSums(is.na(raster::extract(envs, occs))) > 0)
# bg <- dismo::randomPoints(envs[[1]], 10000)
# mod.args <- list(fc = "LQ", rm = 1)
# partitions.kfold <- ENMeval::get.block(occs, bg)
# envs.xy <- raster::rasterToPoints(envs[[1]], spatial = TRUE)
# partitions.kfold$envs.grp <- ENMeval::get.block(occs, envs.xy@coords)$bg.grp
# partitions.split <- ENMeval::get.randomkfold(occs, bg, 2)
# x1 <- nullSDMs(occs, envs, bg, partitions, dismo::maxent, mod.args, 10, "biome")
# x2 <- nullSDMs(occs, envs, bg, partitions.kfold, maxnet::maxnet, mod.args, 5, "kfold", "biome")
# x3 <- nullSDMs(occs, envs, bg, partitions.split, maxnet::maxnet, mod.args, 5, "split", "biome")


# for split evaluation, label training occs "1" and independent evaluation occs "2" in partitions
nullSDMs <- function(occs, envs, bg, partitions,
                     mod.fun, mod.args, no.iter, evalType = c("split", "kfold"),
                     categoricals = NULL, abs.auc.diff = FALSE) {

  if(evalType == "split") {
    # partition number equals 1
    nk <- 1


  }else if(evalType == "kfold") {
    # get total number of partitions (k)
    nk <- length(unique(partitions$occ.grp))
    # define new raster (envs.cv) with values equal to partitions based on partitions$envs.grp
    envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
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

  # get model function name
  mod.fun.name <- as.character(substitute(mod.fun))[3]

  # initialize data frames to collect evaluation stats for each partition
  # per iteration and their averages
  tbl.fields <- c("auc.train", "avg.auc.test", "var.auc.test", "avg.auc.diff",
                  "var.auc.diff", "avg.or.min", "var.or.min", "avg.or.10", "var.or.10")
  k.fields <- c("auc.test", "auc.diff", "or.min", "or.10")
  real.stats <- data.frame(matrix(nrow = 1, ncol = length(tbl.fields)))
  names(real.stats) <- tbl.fields
  null.stats <- data.frame(matrix(nrow = no.iter, ncol = length(tbl.fields)))
  names(null.stats) <- tbl.fields

  ############################## #
  # build real model ####
  ############################## #

  mod.args.real <- buildModArgs(mod.fun.name, mod.args, occs.vals, bg.vals)
  mod.real <- do.call(mod.fun, mod.args.real)
  # calculate training auc
  real.stats$auc.train <- dismo::evaluate(occs.vals, bg.vals, mod.real)@auc
  real.stats$nparam <- no.params(mod.real, mod.fun.name)
  kstats <- data.frame(matrix(nrow = nk, ncol = length(k.fields)))
  names(kstats) <- k.fields

  for(k in 1:nk) {
    if(evalType == "kfold") {
      occs.train.real.k <- occs.vals[partitions$occ.grp != k, ]
      occs.test.real.k <- occs.vals[partitions$occ.grp == k, ]
      bg.train.k <- bg.vals[partitions$bg.grp != k, ]
    }else if(evalType == "split") {
      occs.train.real.k <- occs.vals[partitions$occ.grp != 2,]
      occs.test.real.k <- occs.vals[partitions$occ.grp != 1,]
      bg.train.k <- bg.vals
    }
    mod.args.real.k <- buildModArgs(mod.fun.name, mod.args, occs.train.real.k, bg.train.k)
    mod.real.k <- do.call(mod.fun, mod.args.real.k)
    kstats[k,] <- evalStats(occs.train.real.k, bg.train.k, occs.test.real.k,
                               mod.real.k, abs.auc.diff)
    message(sprintf("Completed real model and evaluation for partition %i.", k))
  }
  # fill in rest of real model statistics
  real.stats$avg.auc.test <- mean(kstats$auc.test)
  real.stats$var.auc.test <- corrected.var(kstats$auc.test, real.stats$nparam)
  real.stats$avg.auc.diff <- mean(kstats$auc.diff)
  real.stats$var.auc.diff <- corrected.var(kstats$auc.diff, real.stats$nparam)
  real.stats$avg.or.min <- mean(kstats$or.min)
  real.stats$var.or.min <- var(kstats$or.min)
  real.stats$avg.or.10 <- mean(kstats$or.10)
  real.stats$var.or.10 <- var(kstats$or.10)

  ############################## #
  # build null models ####
  ############################## #

  # initialize list to record stats for null iteration i
  null.stats.i <- list()

  for(i in 1:no.iter) {
    # initialize data frame for partition statistics for current iteration
    null.stats.i[[i]] <- data.frame(matrix(nrow = nk, ncol = length(k.fields)))
    names(null.stats.i[[i]]) <- k.fields

    # sample random occurrences
    occs.null <- list()
    if(evalType == "kfold") {
      # if kfold evaluation, randomly sample the same number of training
      # occs over each k partition of envs
      for(k in 1:nk) {
        envs.k <- envs.cv
        envs.k[envs.k != k] <- NA
        occs.null.xy <- dismo::randomPoints(envs.k, occs.parts.tbl[k])
        occs.null[[k]] <- as.data.frame(extract(envs, occs.null.xy))
      }
    }else if(evalType == "split") {
      # if split evaluation, randomly sample the same number of training
      # occs over the envs extent only once
      occs.null.xy <- dismo::randomPoints(envs, occs.parts.tbl[1])
      occs.null[[1]] <- as.data.frame(extract(envs, occs.null.xy))
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

    mod.args.i <- buildModArgs(mod.fun.name, mod.args, occs.null.all, bg.vals)
    mod.i <- do.call(mod.fun, mod.args.i)
    null.stats[i,]$auc.train <- dismo::evaluate(occs.null.all, bg.vals, mod.i)@auc
    nparam.i <- no.params(mod.i, mod.fun.name)

    # get training and testing data for current iteration
    for(k in 1:nk) {
      if(evalType == "kfold") {
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
      mod.args.k <- buildModArgs(mod.fun.name, mod.args, occs.null.train, bg.train)
      # build the model
      mod.k <- do.call(mod.fun, mod.args.k)

      null.stats.i[[i]][k,] <- evalStats(occs.null.train, bg.train, occs.test,
                                         mod.k, abs.auc.diff)

      message(sprintf("Completed null model and evaluation for partition %i of iteration %i.", k, i))
    }
    # average partition statistics
    null.stats[i,]$avg.auc.test <- mean(null.stats.i[[i]]$auc.test)
    null.stats[i,]$var.auc.test <- corrected.var(null.stats.i[[i]]$auc.test, nparam.i)
    null.stats[i,]$avg.auc.diff <- mean(null.stats.i[[i]]$auc.diff)
    null.stats[i,]$var.auc.diff <- corrected.var(null.stats.i[[i]]$auc.diff, nparam.i)
    null.stats[i,]$avg.or.min <- mean(null.stats.i[[i]]$or.min)
    null.stats[i,]$var.or.min <- var(null.stats.i[[i]]$or.min)
    null.stats[i,]$avg.or.10 <- mean(null.stats.i[[i]]$or.10)
    null.stats[i,]$var.or.10 <- var(null.stats.i[[i]]$or.10)
  }

  out <- list(real.stats = real.stats, null.stats = null.stats, null.stats.i = null.stats.i)
  return(out)
}

#dismo::maxent(x, p, args = c(args.i, userArgs),factors = categoricals)
