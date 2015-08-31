

read_stats <- function() {
  stats1 <<- read.table("seeds1.stats.txt", header = TRUE)
  stats2 <<- read.table("seeds2.stats.txt", header = TRUE)
}

plot.sen.by.len <- function(org, T) {
  T <- T[T$org == org,]

  plot(c(), c(), xlim = c(0,max(T$l, na.rm=TRUE)), ylim = c(0, max(c(T$tpr,T$QuCoverage))),
       xlab = "Seed length", ylab = "sensitivity", main = org)
  abline(h=T[T$tool == "RepeatScout","tpr"], col = 'red')
  abline(h=T[T$tool == "RepeatScout","QuCoverage"], col = 'red', lty = 'dashed')
  
  T2 = T[T$tool == 'phRAIDER',]
  lines(T2$l, T2$tpr, col='blue')
  lines(T2$l, T2$QuCoverage, col = 'blue', lty = 'dashed')
}

plot.sen.by.ratio <- function(org, T) {
  T <- T[T$org == org,]
  
  plot(c(), c(), xlim = c(0,max(T$w.l, na.rm=TRUE)), ylim = c(0, max(c(T$tpr,T$QuCoverage), na.rm=TRUE)),
       xlab = "Seed length", ylab = "sensitivity", main = org)
  abline(h=T[T$tool == "RepeatScout","tpr"], col = 'red')
  abline(h=T[T$tool == "RepeatScout","QuCoverage"], col = 'red', lty = 'dashed')
  
  T2 = T[T$tool == 'phRAIDER',]
  points(T2$w.l, T2$tpr, col='blue')
  points(T2$w.l, T2$QuCoverage, col = 'green')
}

fixed_plots <- function() {
  quartz()
  plot.sen.by.len("hg38.chr22", stats)
  quartz()
  plot.sen.by.len("ce10.chrV", stats)
}