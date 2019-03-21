#' @export

evalStats <- function(occs.train, bg.train, occs.test, mod, abs.auc.diff) {
  # calculate auc on training and testing data
  auc.train <- dismo::evaluate(occs.train, bg.train, mod)@auc
  auc.test <- dismo::evaluate(occs.test, bg.train, mod)@auc
  # calculate auc diff
  auc.diff <- auc.train - auc.test
  if(abs.auc.diff == TRUE) auc.diff <- abs(auc.diff)
  # get model predictions for training and testing data
  pred.train <- dismo::predict(mod, occs.train)
  pred.test <- dismo::predict(mod, occs.test)
  # get 10 percentile predicted value
  occs.train.n <- nrow(occs.train)
  if(occs.train.n < 10) {
    pct10.train <- floor(occs.train.n * 0.1)
  } else {
    pct10.train <- ceiling(occs.train.n * 0.1)
  }
  pct10.train.thr <- sort(pred.train)[pct10.train]
  or10.test <- mean(pred.test < pct10.train.thr)
  min.train.thr <- min(pred.train)
  orMin.test <- mean(pred.test < min.train.thr)

  stats <- c(auc.test, auc.diff, orMin.test, or10.test)

  return(stats)
}
