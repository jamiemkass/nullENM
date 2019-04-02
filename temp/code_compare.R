# corentin code begins with c
# nullENM code begins with n

########################################## #
# run Corentin's code ####
########################################## #
set.seed(48)
setwd("/Users/musasabi/Documents/research/nullSDM_testing")
bg.files <- list.files(path="data", pattern="background", full.names=TRUE)
bgs <- lapply(bg.files, read.csv)
AR.bg <- bgs[[1]][,4:27]
ES.bg <- bgs[[3]][,4:27]
FLs.bg <- bgs[[5]][,4:27]
NY.bg <- bgs[[8]][,4:27]
ARESNYFLs.bg <- bgs[[2]][,4:27]
TEST.bg <- bgs[[6]]
global.bg <- bgs[[7]][,4:27]

bgs <- list(ARESNYFLs.bg, AR.bg, ES.bg, NY.bg, FLs.bg)

train.files <- list.files(path="data", pattern="train", full.names=TRUE)
trains <- lapply(train.files, read.csv)
AR.train <- trains[[1]][,4:27]
ARESNYFLs.train <- trains[[2]][,4:27]
ARESNYFLs.test <- read.csv('data/MOPAm_noout_ARESNYFLs_test.csv')[,4:27]
ARESNYFLs <- rbind(ARESNYFLs.train, ARESNYFLs.test)
ES.train <- trains[[3]][,4:27]
FLs.train <- trains[[4]][,4:27]
NY.train <- trains[[5]][,4:27]

trains <- list(ARESNYFLs, AR.train, ES.train, NY.train, FLs.train)
names(trains) <- c("ARESNYFLs", "AR", "ES", "NY", "FLs")

test <- read.csv('data/MOPAmnout_global_test.csv')
test <- unique(test[,4:27])

source("/Users/musasabi/Documents/github/nullSDM/temp/corentin_nullSDMs.R")

#run models
library(rJava) # required to run Maxent within dismo
options(java.parameters = "-Xmx1g" ) # optional, to set the memory allocated to java, hence maxent has to be done before loading dismo
library(dismo) # should also load automatically the required packages sp and raster

# default settings
c.default <- list()
for(i in 1:length(trains)) {
  x <- mxt.nulltest(trains[[i]], bgs[[i]], n = 100, proj=bgs[[i]], c.proj=TRUE, args="noaddsamplestobackground", test=test)
  c.default[[names(trains)[i]]] <- x
}
saveRDS(c.default, "results/c.default.rds")
c.default <- readRDS("results/c.default.rds")

# simplest settings
args=c("noaddsamplestobackground","noautofeature","noproduct","nothreshold","nohinge","noquadratic","betamultiplier=4")
c.L4 <- list()
for(i in 1:length(trains)) {
  x <- mxt.nulltest(trains[[i]], bgs[[i]], n = 100, proj=global.bg, c.proj=TRUE, args=args, test=test)
  c.L4[[names(trains)[i]]] <- x
}
saveRDS(c.L4, "results/c.L4.rds")

# most complex settings
args=c("noaddsamplestobackground","noautofeature","betamultiplier=0.25")
c.LQPTH025 <- list()
for(i in 1:length(trains)) {
  x <- mxt.nulltest(trains[[i]], bgs[[i]], n = 100, proj=global.bg, c.proj=TRUE, args=args, test=test)
  c.LQPTH025[[names(trains)[i]]] <- x
}
saveRDS(c.LQPTH025, "results/c.LQPTH025.rds")

########################################## #
# run nullENM code ####
########################################## #

# default settings
n.default <- list()
for(i in 1:length(trains)) {
  occs <- trains[[i]]
  bg <- bgs[[i]]
  x <- nullSDMs(occs = occs, bg = bg, occs.indTest = test,
                mod.name = "maxent.jar", mod.args = list(fc = "LQHPT", rm = 1),
                no.iter = 100, eval.type = "split")
  n.default[[names(trains)[i]]] <- x
}
saveRDS(n.default, "results/n.default.rds")

# simplest settings
n.L4 <- list()
for(i in 1:length(trains)) {
  occs <- rbind(trains[[i]], test)
  occs.grp <- c(rep(1, nrow(trains[[i]])), rep(2, nrow(test)))
  bg.grp <- rep(0, nrow(bgs[[i]]))
  x <- nullSDMs(occs = occs, bg = bgs[[i]], occs.grp = occs.grp, bg.grp = bg.grp,
                mod.name = "maxent.jar", mod.args = list(fc = "L", rm = 4),
                no.iter = 100, eval.type = "split")
  n.L4[[names(trains)[i]]] <- x
}
saveRDS(n.L4, "results/n.L4.rds")

# most complex settings
n.LQPTH025 <- list()
for(i in 1:length(trains)) {
  occs <- rbind(trains[[i]], test)
  occs.grp <- c(rep(1, nrow(trains[[i]])), rep(2, nrow(test)))
  bg.grp <- rep(0, nrow(bgs[[i]]))
  x <- nullSDMs(occs = occs, bg = bgs[[i]], occs.grp = occs.grp, bg.grp = bg.grp,
                mod.name = "maxent.jar", mod.args = list(fc = "LQHPT", rm = 0.25),
                no.iter = 100, eval.type = "split")
  n.LQPTH025[[names(trains)[i]]] <- x
}
saveRDS(n.LQPTH025, "results/n.LQPTH025.rds")

