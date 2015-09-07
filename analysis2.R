read.file <- function(file, pal, test_num) {
  T <- read.table(file, header = TRUE)
  T$pal = pal
  T$test_num = test_num
  return(T)
}

read.stats <- function() {
  stats1 <<- read.file("seeds1.stats.txt", FALSE, 1) # Vary length
  stats2 <<- read.file("seeds2.stats.txt", FALSE, 2)
  stats3 <<- read.file("seeds3.stats.txt", FALSE, 3)
  stats4 <<- read.file("seeds4.stats.txt", FALSE, 4)
  stats5 <<- read.file("seeds5.stats.txt", TRUE,  5)
  stats6 <<- read.file("seeds6.stats.txt", FALSE, 6)
  stats7 <<- read.file("seeds7.stats.txt", TRUE,  7)
  stats8 <<- read.file("seeds8.stats.txt", FALSE, 8)
  stats9 <<- read.file("seeds9.stats.txt", TRUE,  9)
  stats10 <<- read.file("seeds10.stats.txt", FALSE, 10)
  stats11 <<- read.file("seeds11.stats.txt", TRUE,  11)  
  stats.all <<-rbind(stats1, stats2, stats3, stats4, stats5, stats6, stats7, stats8, stats9, stats10, stats11)
}

add.legend <- function(position) {
  legend(position, c("RM", "PRA"), col = c('blue', 'green'), 
                     lty = c('solid', 'solid', 'solid', 'solid'))
}

plot.gen <- function(T, org, x_col, y_col, color, xlab = NA, ylab = NA, main = NA, 
                     xlim = NA, ylim = NA, legend_position = "topleft", new_plot = TRUE,
                     reverse_x = FALSE, reverse_y = FALSE) {
  T <- droplevels(T[T$org==org,])

  if (reverse_x) {
    xlim = rev(xlim)
  }
  if (reverse_y) {
    ylim = rev(ylim)
  }
  
  if (new_plot) {
    plot(c(), c(), xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  }
  abline(h = subset(T, tool == "RepeatScout", select = y_col)[[1]], col = color)
  points(subset(T, tool=phRAIDER, select=c(x_col, y_col)), col = color)
}
  
  
plot.ROC <- function(T, org, tpr, tnr, pivot, color) {
  T <- T[T$org == org & T$tool == 'phRAIDER' & !is.na(T[tpr]) & !is.na(T[tnr]), ]
  X <- c()
  Y <- c()
  for (v in unique(T[[pivot]])) {
    T2 = T[T[pivot] == v,]
    Y <- c(Y, mean(T2[[tpr]], na.rm=TRUE))
    X <- c(X, 1 - mean(T2[[tnr]]))
  }
  plot(X, Y, xlim = c(0,1), ylim = c(0,1), col=color, main = "ROC curve", xlab = "FNR", ylab = "TPR")
}


plot.test <- function(T, org, x_axis = "w.l", x_label = "w/l") {
  par(mfrow=c(1,2))
  # Sensitity v. weight/length
  xlim = range(T[T$org==org, x_axis], na.rm=TRUE)
  ylim = range(T[T$org==org, c('tpr', 'QuCoverage')], na.rm=TRUE)
  plot.gen(T, org, x_axis, "tpr", "blue", xlab = x_label, ylab = "", xlim =xlim,
           ylim = ylim, main = "Sensitivity", new_plot = TRUE, reverse_x = FALSE)
  plot.gen(T, org, x_axis, "QuCoverage", "red", new_plot = FALSE, reverse_x = FALSE)
  #add.legend("bottomleft")
  
  # Specificity v. weight/length
  xlim = range(T[T$org==org, x_axis], na.rm=TRUE)
  ylim = range(T[T$org==org, c('tnr', 'ConCoverage')], na.rm=TRUE)
  plot.gen(T, org, x_axis, "tnr", "blue", xlab =  x_label, ylab = "", xlim =xlim,
           ylim = ylim, main = "Specificity", new_plot = TRUE, reverse_x = FALSE)
  plot.gen(T, org, x_axis, "ConCoverage", "red", new_plot = FALSE, reverse_x = FALSE) 
  
#   # Runtime v. weight/length
#   xlim = range(T[T$org==org, "w.l"], na.rm=TRUE)
#   ylim = range(T[T$org==org, c('ToolCpuTime', 'ConCoverage')], na.rm=TRUE)
#   plot.gen(T, org, "w.l", "ToolCpuTime", "blue", xlab = "w/l", ylab = "runtime", xlim = xlim,
#            ylim = ylim, main = "Runtime", new_plot = TRUE, reverse_x = TRUE)
#   
#   # ROC
#   plot.ROC(T, org, "tpr", "tnr", "w.l", "blue")
#   #plot.ROC(T, org, "QuCoverage", "ConCoverage", "w.l", "green")

  
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