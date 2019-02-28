plot.nullSDMs <- function(results, stat) {
  null.tbl.name <- ifelse(stat == "auc.train", "auc.train", paste0("mean.", stat))
  null.stats <- round(results@null.stats[[null.tbl.name]], 3)
  real.stat <- round(results@all.stats["real.mean", stat], 3)
  stat.max <- ifelse(real.stat < max(null.stats), max(null.stats), real.stat)
  hist(null.stats, xlim=c(min(null.stats), stat.max + stat.max/10), main = "", xlab="AUC test")
  abline(v=quantile(null.stats, 0.95), col='blue', lty=3)
  abline(v=quantile(null.stats, 0.99), col='blue', lty=2)
  abline(v=real.stat, col='red')
}
