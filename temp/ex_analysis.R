# bv <- spocc::occ('Bradypus variegatus', 'gbif', limit=300, has_coords=TRUE)
# occs <- as.data.frame(bv$gbif$data$Bradypus_variegatus[,2:3])
# occs <- occs[!duplicated(occs),]
# envs <- raster::stack(list.files(path=paste(system.file(package='dismo'), '/ex', sep=''), pattern='grd', full.names=TRUE))
# envs <- raster::mask(envs, envs$biome)
# which(rowSums(is.na(raster::extract(envs, occs))) > 0)
# bg <- dismo::randomPoints(envs[[1]], 10000)
# mod.args <- list(fc = "LQ", rm = 1)
# p.block <- ENMeval::get.block(occs, bg)
# envs.xy <- raster::rasterToPoints(envs[[1]], spatial = TRUE)
# p.block$envs.grp <- ENMeval::get.block(occs, envs.xy@coords)$bg.grp
# p.rand <- ENMeval::get.randomkfold(occs, bg, 3)
# p.split <- ENMeval::get.randomkfold(occs, bg, 2)
# x1 <- nullSDMs(occs, envs, bg, p.rand$occ.grp, p.rand$bg.grp, "maxnet", mod.args, 5, "kfold", "biome")
# x2 <- nullSDMs(occs, envs, bg, p.block$occ.grp, p.block$bg.grp, envs.grp = p.block$envs.grp, "maxnet", mod.args, 5, "kspatial", "biome")
# x3 <- nullSDMs(occs, envs, bg, p.split$occ.grp, p.split$bg.grp, "maxnet", mod.args, 5, "split", "biome")

bg.files=list.files(path="/Users/jamie/Documents/nullSDM_testing/corentin_data",pattern="background",full.names=TRUE)
bgs=lapply(bg.files,read.csv)
AR.bg=bgs[[1]][,4:27]
ES.bg=bgs[[3]][,4:27]
FLs.bg=bgs[[5]][,4:27]
NY.bg=bgs[[8]][,4:27]
ARESNYFLs.bg=bgs[[2]][,4:27]
TEST.bg=bgs[[6]]
global.bg=bgs[[7]][,4:27]

BGl=list(ARESNYFLs.bg,AR.bg,ES.bg,NY.bg,FLs.bg)

rm(bg.files,bgs)

train.files=list.files(path="/Users/jamie/Documents/nullSDM_testing/corentin_data",pattern="train",full.names=TRUE)
trains=lapply(train.files,read.csv)
AR.train=trains[[1]][,4:27]
ARESNYFLs.train=trains[[2]][,4:27]
ARESNYFLs.test=read.csv('/Users/jamie/Documents/nullSDM_testing/corentin_data/MOPAm_noout_ARESNYFLs_test.csv')[,4:27]
ARESNYFLs=rbind(ARESNYFLs.train,ARESNYFLs.test)
ES.train=trains[[3]][,4:27]
FLs.train=trains[[4]][,4:27]
NY.train=trains[[5]][,4:27]

TRAINl=list(ARESNYFLs,AR.train,ES.train,NY.train,FLs.train)

rm(train.files,trains,ARESNYFLs.train,ARESNYFLs.test)

test=read.csv('/Users/jamie/Documents/nullSDM_testing/corentin_data/MOPAmnout_global_test.csv')
test=unique(test[,4:27])

SAVEl=lapply(list('ARESNYFLs','AR','ES','NY','FLs'),function(x) {paste('/Users/jamie/Documents/nullSDM_testing/corentin_data/null_runs/MOPApred_',x,"_default.csv",sep="")})

corentin_runs <- list()
cl <- parallel::makeCluster(4)
doSNOW::registerDoSNOW(cl)
library(foreach)

mod.args.lst <- list(list(fc = "LQHPT", rm = 1),
                     list(fc = "L", rm = 4),
                     list(fc = "LQHPT", rm = 0.25))

corentin_runs <- foreach(i = 1:length(TRAINl), .packages = c("dismo", "raster", "rJava", "nullSDM")) %dopar% {
  occs <- rbind(TRAINl[[i]], test)
  occs.grp <- c(rep(1, nrow(TRAINl[[i]])), rep(2, nrow(test)))
  bg.grp <- rep(0, nrow(BGl[[i]]))
  lst <- list()
  for(j in 1:3) {
    lst[[j]] <- nullSDMs(occs=occs, bg=BGl[[i]],
                         occs.grp = occs.grp, bg.grp = bg.grp,
                         mod.name = "maxent.jar", mod.args = mod.args.lst[[j]], no.iter = 5,
                         eval.type = "split")
  }
  lst
}
parallel::stopCluster(cl)
