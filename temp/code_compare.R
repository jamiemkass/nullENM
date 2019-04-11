# corentin code begins with c
# nullENM code begins with n

########################################## #
# run Corentin's code ####
########################################## #
set.seed(48)
setwd("/Users/kass/Documents/nullENM_testing")
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

source("/Users/kass/Documents/github/nullENM/temp/corentin_nullSDMs.R")

#run models
library(rJava) # required to run Maxent within dismo
options(java.parameters = "-Xmx1g" ) # optional, to set the memory allocated to java, hence maxent has to be done before loading dismo
library(dismo) # should also load automatically the required packages sp and raster

# default settings
c.default <- list()
for(i in 1:length(trains)) {
  x <- mxt.nulltest(trains[[i]], bgs[[i]], n = 1000, proj=bgs[[i]], c.proj=TRUE, args="noaddsamplestobackground", test=test)
  c.default[[names(trains)[i]]] <- x
}
saveRDS(c.default, "results/c.default.rds")
c.default <- readRDS("results/c.default.rds")

# simplest settings
args=c("noaddsamplestobackground","noautofeature","noproduct","nothreshold","nohinge","noquadratic","betamultiplier=4")
c.L4 <- list()
for(i in 1:length(trains)) {
  x <- mxt.nulltest(trains[[i]], bgs[[i]], n = 1000, proj=bgs[[i]], c.proj=TRUE, args=args, test=test)
  c.L4[[names(trains)[i]]] <- x
}
saveRDS(c.L4, "results/c.L4.rds")
c.L4 <- readRDS("results/c.L4.rds")

# most complex settings
args=c("noaddsamplestobackground","noautofeature","betamultiplier=0.25")
c.LQPTH025 <- list()
for(i in 1:length(trains)) {
  x <- mxt.nulltest(trains[[i]], bgs[[i]], n = 1000, proj=bgs[[i]], c.proj=TRUE, args=args, test=test)
  c.LQPTH025[[names(trains)[i]]] <- x
}
saveRDS(c.LQPTH025, "results/c.LQPTH025.rds")
c.LQPTH025 <- readRDS("results/c.LQPTH025.rds")

########################################## #
# run nullENM code ####
########################################## #

# default settings
n.default <- list()
for(i in 1:length(trains)) {
  occs <- trains[[i]]
  bg <- bgs[[i]]
  x <- nullENMs(occs = occs, bg = bg, occs.indTest = test,
                mod.name = "maxent.jar", mod.args = list(fc = "LQHPT", rm = 1),
                no.iter = 1000, eval.type = "split")
  n.default[[names(trains)[i]]] <- x
}
saveRDS(n.default, "results/n.default.rds")
n.default <- readRDS("results/n.default.rds")

# simplest settings
n.L4 <- list()
for(i in 1:length(trains)) {
  occs <- trains[[i]]
  bg <- bgs[[i]]
  x <- nullENMs(occs = occs, bg = bg, occs.indTest = test,
                mod.name = "maxent.jar", mod.args = list(fc = "L", rm = 4),
                no.iter = 1000, eval.type = "split")
  n.L4[[names(trains)[i]]] <- x
}
saveRDS(n.L4, "results/n.L4.rds")
n.L4 <- readRDS("results/n.L4.rds")

# most complex settings
n.LQPTH025 <- list()
for(i in 1:length(trains)) {
  occs <- trains[[i]]
  bg <- bgs[[i]]
  x <- nullENMs(occs = occs, bg = bg, occs.indTest = test,
                mod.name = "maxent.jar", mod.args = list(fc = "LQHPT", rm = 0.25),
                no.iter = 1000, eval.type = "split")
  n.LQPTH025[[names(trains)[i]]] <- x
}
saveRDS(n.LQPTH025, "results/n.LQPTH025.rds")
n.LQPTH025 <- readRDS("results/n.LQPTH025.rds")

########################################## #
# run nullENM code ####
########################################## #
library(dplyr)
library(tidyr)
sum.stat.names <- rownames(n.default$ARESNYFLs@all.stats)
c.default <- lapply(c.default, function(x) x@summary %>% select(auc.train = AUCtrain, auc.test = AUCtest, auc.diff = AUCdiff, or.10 = OR) %>% mutate(code = "C", run = "default"))
c.L4 <- lapply(c.L4, function(x) x@summary %>% select(auc.train = AUCtrain, auc.test = AUCtest, auc.diff = AUCdiff, or.10 = OR) %>% mutate(code = "C", run = "L4"))
c.LQPTH025 <- lapply(c.LQPTH025, function(x) x@summary %>% select(auc.train = AUCtrain, auc.test = AUCtest, auc.diff = AUCdiff, or.10 = OR) %>% mutate(code = "C", run = "LQPTH025"))
c.all <- rbind(do.call(rbind, c.default), do.call(rbind, c.L4), do.call(rbind, c.LQPTH025))
c.all$region <- rep(names(trains), 3, each = 5)
c.all$sum.stat <- rep(sum.stat.names[-2], 15)
row.names(c.all) <- NULL

n.default <- lapply(n.default, function(x) x@all.stats %>% select(auc.train, auc.test, auc.diff, or.10) %>% mutate(code = "N", run = "default"))
n.L4 <- lapply(n.L4, function(x) x@all.stats %>% select(auc.train, auc.test, auc.diff, or.10) %>% mutate(code = "N", run = "L4"))
n.LQPTH025 <- lapply(n.LQPTH025, function(x) x@all.stats %>% select(auc.train, auc.test, auc.diff, or.10) %>% mutate(code = "N", run = "LQPTH025"))
n.all <- rbind(do.call(rbind, n.default), do.call(rbind, n.L4), do.call(rbind, n.LQPTH025))
n.all$region <- rep(names(trains), 3, each = 6)
n.all$sum.stat <- rep(sum.stat.names, 15)
row.names(n.all) <- NULL
n.all <- na.omit(n.all)

df <- data.frame(as.matrix(c.all[1:4]) - as.matrix(n.all[1:4]), c.all[6:8])
df <- df %>% gather(eval.stat, value, auc.train:or.10)

library(ggplot2)

ggplot(df, aes(x = eval.stat, y = value, color = region)) +
  facet_grid(rows = vars(sum.stat), scales = "free_y") + geom_point()
