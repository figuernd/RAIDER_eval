read_stats <- function(file, pal, test_num) {
  T <- read.table(file, header = TRUE)
  T$pal = pal
  T$test_num = test_num
  return(T)
}

read.stats <- function() {
  stats1 <<- read_stats("seeds1.stats.txt", FALSE, 1)
  stats2 <<- read_stats("seeds2.stats.txt", FALSE, 2)
  stats3 <<- read_stats("seeds3.stats.txt", FALSE, 3)
  stats4 <<- read_stats("seeds4.stats.txt", FALSE, 4)
  stats5 <<- read_stats("seeds5.stats.txt", TRUE,  5)
  stats6 <<- read_stats("seeds6.stats.txt", FALSE, 6)
  stats7 <<- read_stats("seeds7.stats.txt", TRUE,  7)
}

plot.sen.by.len <- function(org, T) {
  T <- T[T$org == org,]

  plot(c(), c(), xlim = range(T$l, na.rm=TRUE), ylim = range(c(T$tpr,T$QuCoverage)),
       xlab = "Seed length", ylab = "sensitivity", main = org)
  legend("topleft", c("RptScout (RM)", "RptScout (PRA)", "RAIDER (RM)", "RAIDER (PRM)"), col = c('red', 'red', 'blue', 'green'), lty = c('solid', 'dashed', 'solid', 'solid'))
  abline(h=T[T$tool == "RepeatScout","tpr"], col = 'red')
  abline(h=T[T$tool == "RepeatScout","QuCoverage"], col = 'red', lty = 'dashed')
  
  T2 = T[T$tool == 'phRAIDER',]
  points(T2$l, T2$tpr, col='blue')
  points(T2$l, T2$QuCoverage, col = 'green', lty = 'dashed')
}

plot.sen.by.ratio <- function(org, T) {
  T <- T[T$org == org,]
  print(T[1,"test_num"])
  plot(c(), c(), xlim = range(T$w.l, na.rm=TRUE), ylim = range(c(T$tpr,T$QuCoverage), na.rm=TRUE),
       xlab = "W/L ratio", ylab = "sensitivity", main = sprintf("test%d", T[1,"test_num"]))
  abline(h=T[T$tool == "RepeatScout","tpr"], col = 'red')
  abline(h=T[T$tool == "RepeatScout","QuCoverage"], col = 'red', lty = 'dashed')
  
  T2 = T[T$tool == 'phRAIDER',]
  points(T2$w.l, T2$tpr, col='blue')
  points(T2$w.l, T2$QuCoverage, col = 'green')
}

fixed_plots <- function() {
  par(mfrow=c(3,3))
  plot.sen.by.len("ce10.chrV", stats1)
  plot.sen.by.ratio("ce10.chrV", stats2)
  plot.sen.by.ratio("ce10.chrV", stats3)
  plot.sen.by.ratio("ce10.chrV", stats4)
  plot.sen.by.ratio("ce10.chrV", stats5)
  plot.sen.by.ratio("ce10.chrV", stats6)
  plot.sen.by.ratio("ce10.chrV", stats7)
}