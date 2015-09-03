read.stats <- function() {
  stats1 <<- read.table("seeds1.stats.txt", header = TRUE)
  stats1$seed = NA
  stats1$note = NA
  stats2 <<- read.table("seeds2.stats.txt", header = TRUE)
  stats2$seed = NA
  stats2$note = NA
  stats2 <<- stats2[stats2$w.l > 0, ]
  stats3 <<- read.table("seeds3.stats.txt", header = TRUE)
  stats3$note = NA
  stats4 <<- read.table("seeds4.stats.txt", header = TRUE)
  stats4$note = NA
  stats5 <<- read.table("seeds5.stats.txt", header = TRUE)
  stats4$note = "pal"
}

plot.sen.by.len <- function(org, T) {
  T <- T[T$org == org,]

  plot(c(), c(), xlim = range(T$l, na.rm=TRUE), ylim = range(c(T$tpr,T$QuCoverage)),
       xlab = "Seed length", ylab = "sensitivity", main = org)
  abline(h=T[T$tool == "RepeatScout","tpr"], col = 'red')
  abline(h=T[T$tool == "RepeatScout","QuCoverage"], col = 'red', lty = 'dashed')
  
  T2 = T[T$tool == 'phRAIDER',]
  points(T2$l, T2$tpr, col='blue')
  points(T2$l, T2$QuCoverage, col = 'green', lty = 'dashed')
}

plot.sen.by.ratio <- function(org, T) {
  T <- T[T$org == org,]
  
  plot(c(), c(), xlim = range(T$w.l, na.rm=TRUE), ylim = range(c(T$tpr,T$QuCoverage), na.rm=TRUE),
       xlab = "W/L ratio", ylab = "sensitivity", main = org)
  abline(h=T[T$tool == "RepeatScout","tpr"], col = 'red')
  abline(h=T[T$tool == "RepeatScout","QuCoverage"], col = 'red', lty = 'dashed')
  
  T2 = T[T$tool == 'phRAIDER',]
  points(T2$w.l, T2$tpr, col='blue')
  points(T2$w.l, T2$QuCoverage, col = 'green')
}

fixed_plots <- function() {
  par(mfrow=c(3,2))
  plot.sen.by.len("ce10.chrV", stats1)
  plot.sen.by.ratio("ce10.chrV", stats2)
  plot.sen.by.len("ce10.chrV", stats3)
  plot.sen.by.ratio("ce10.chrV", stats3)
  plot.sen.by.len("ce10.chrV", stats4)
  plot.sen.by.ratio("ce10.chrV", stats4)
  plot.sen.by.len("ce10.chrV", stats5)
  plot.sen.by.ratio("ce10.chrV", stats5)
}