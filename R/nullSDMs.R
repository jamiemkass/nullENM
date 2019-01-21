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
# partitions.kfold$occs <- partitions.kfold$occ.grp
# partitions.kfold$bg <- partitions.kfold$bg.grp
# envs.xy <- raster::rasterToPoints(envs[[1]], spatial = TRUE)
# partitions.kfold$envs <- ENMeval::get.block(occs, envs.xy@coords)$bg.grp
# partitions.split <- ENMeval::get.randomkfold(occs, bg, 2)
# x1 <- nullSDMs(occs, envs, bg, partitions, dismo::maxent, mod.args, 10, "biome")
# x2 <- nullSDMs(occs, envs, bg, partitions.kfold, maxnet::maxnet, mod.args, 5, "kfold", "biome")

nullSDMs <- function(occs, envs, bg, partitions,
                     mod.fun, mod.args, no.iter, eval = c("split", "kfold"),
                     categoricals = NULL, abs.auc.diff = FALSE) {

  if(eval == "split") {
    # partition number equals 1
    nk <- 1


  }

  if(eval == "kfold") {
    # get total number of partitions (k)
    nk <- length(unique(partitions$occs))
    # define new raster (envs.cv) with values equal to partitions based on partitions$envs
    envs.xy <- rasterToPoints(envs[[1]], spatial = TRUE)
    envs.cellNo <- raster::extract(envs[[1]], envs.xy, cellnumbers = TRUE)
    envs.cv <- envs[[1]]
    envs.cv[envs.cellNo[,1]] <- partitions$envs
  }

  # get number of occurrence points by partition
  occs.parts.tbl <- table(partitions$occs)

  # get environmental values for occs and bg
  occs.vals <- as.data.frame(raster::extract(envs, occs))
  bg.vals <- as.data.frame(raster::extract(envs, bg))

  # get model function name
  mod.fun.name <- as.character(substitute(mod.fun))[3]

  # initialize data frames to collect evaluation stats for each partition
  # per iteration and their averages
  tbl.fields <- c("auc.train", "avg.auc.test", "std.auc.test", "avg.auc.diff",
                  "std.auc.diff", "avg.or.min", "std.or.min", "avg.or.10", "std.or.10")
  k.fields <- c("auc.test", "auc.diff", "or.min", "or.10")
  real.stats <- data.frame(matrix(nrow = 1, ncol = length(tbl.fields)))
  names(real.stats) = tbl.fields
  null.stats <- data.frame(matrix(nrow = no.iter, ncol = length(tbl.fields)))
  names(null.stats) = tbl.fields



  ############################## #
  # build real model ####
  ############################## #

  mod.args.real <- buildModArgs(mod.fun.name, mod.args, occs.vals, bg.vals)
  mod.real <- do.call(mod.fun, mod.args.real)
  # calculate training auc
  real.stats$auc.train <- dismo::evaluate(occs.vals, bg.vals, mod.real)@auc
  kstats <- data.frame(matrix(nrow = nk, ncol = length(k.fields)))
  names(kstats) <- k.fields

  for(k in 1:nk) {
    occs.train.real.k <- occs.vals[partitions$occs != k, ]
    occs.test.real.k <- occs.vals[partitions$occs == k, ]
    bg.train.k <- bg.vals[partitions$bg != k, ]

    mod.args.real.k <- buildModArgs(mod.fun.name, mod.args, occs.train.real.k, bg.train.k)
    mod.real.k <- do.call(mod.fun, mod.args.real.k)
    kstats[k,] <- evalStats(occs.train.real.k, bg.train.k, occs.test.real.k,
                               mod.real.k, abs.auc.diff)
    message(sprintf("Completed real model and evaluation for partition %i.", k))
  }
  # fill in rest of real model statistics
  real.stats$avg.auc.test <- mean(kstats$auc.test)
  real.stats$std.auc.test <- sd(kstats$auc.test)
  real.stats$avg.auc.diff <- mean(kstats$auc.diff)
  real.stats$std.auc.diff <- sd(kstats$auc.diff)
  real.stats$avg.or.min <- mean(kstats$or.min)
  real.stats$std.or.min <- sd(kstats$or.min)
  real.stats$avg.or.10 <- mean(kstats$or.10)
  real.stats$std.or.10 <- sd(kstats$or.10)



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
    for(k in 1:nk) {
      envs.k <- envs.cv
      envs.k[envs.k != k] <- NA
      occs.null.xy <- dismo::randomPoints(envs.k, occs.parts.tbl[k])
      occs.null[[k]] <- as.data.frame(extract(envs, occs.null.xy))
    }
    occs.null.all <- do.call(rbind, occs.null)

    mod.args.i <- buildModArgs(mod.fun.name, mod.args, occs.null.all, bg.vals)
    mod.i <- do.call(mod.fun, mod.args.i)
    auc.train.i <- dismo::evaluate(occs.null.all, bg.vals, mod.i)@auc

    for(k in 1:nk) {
      # get training and testing data for current iteration
      # (non-k training and k testing partitions)
      occs.null.train <- do.call(rbind, occs.null[-k])
      bg.train <- bg.vals[partitions$bg != k, ]
      occs.test <- occs.vals[partitions$occs == k, ]

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

      null.stats.i[[i]][k,] <- evalStats(occs.null.train, bg.train, occs.test,
                                         mod.k, abs.auc.diff)

      message(sprintf("Completed null model and evaluation for partition %i of iteration %i.", k, i))
    }
    # average partition statistics
    means <- apply(null.stats.i[[i]], 2, mean)
    stds <- apply(null.stats.i[[i]], 2, sd)
    # alternate the above values for data frame
    stats <- c(rbind(means, stds))
    null.stats[i,] <- c(auc.train.i, stats)
  }

  out <- list(real.stats = real.stats, null.stats = null.stats, null.stats.i = null.stats.i)
  return(out)

}

#dismo::maxent(x, p, args = c(args.i, userArgs),factors = categoricals)
