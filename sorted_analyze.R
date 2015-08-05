#!/software/R/3.0.2-gcc/bin/Rscript

# sorted_analyze.R : helper R script for formatting stats data from RAIDER_eval.
#   provides functions for sorting data, as well as adding extra/removing excess
#   information to/from the data. also includes function for outputting data to
#   file/stdout.
# by Carly Schaeffer

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(gdata))

# Reads in data frame from stats file
getDf <- function(fname) {
  dat <- read.table(fname,header=TRUE,comment.char="",fill=TRUE)
  dat <- dat[rowSums(is.na(dat)) < 4, ]
}

sumOverMatches <- function(matches){
  sum(as.numeric(unlist(matches)))
}

getMatches <- function(expr, seed){
  regmatches(seed, gregexpr(expr, seed, perl=TRUE))
}

# Determines number of times val is to be seen in expanded seed (accounting for formatting)
numberOf <- function(seed, val){
  match.string <- paste("(?<=", val, "\\^\\{)[[:digit:]]+(?=\\})", sep="")
  numOnes <- sumOverMatches(getMatches(match.string, as.character(seed)))
  match.string <- paste(val, "(?!\\^)", sep="")
  numOnes + length(unlist(getMatches(match.string, seed)))
}

# Attempts to expand seed that is given in 1^{}, 0^{} format, representing it as a
# regular string. i.e. 1^{3}0^{2}1 ---> 111001
toSeedString <- function(seed){
  seed.string = ""
  match.string <- "(?<=\\^\\{)[[:digit:]]+(?=\\})"
  match.indices <- gregexpr(match.string,seed,perl=TRUE)[[1]]
  match.values <- regmatches(seed,gregexpr(match.string,seed,perl=TRUE))[[1]]
  j = 1
  if (length(match.indices) == 1 && match.indices[1]==-1){
    seed
  }
  else{
    for( i in 1:length(match.indices)){
      next.index <- match.indices[i]
      next.val <- match.values[i]
      seed.val.index = next.index - 3
      if(seed.val.index > j + 1){
        seed.string <- paste(seed.string, substring(seed, j, seed.val.index - 1), sep="")
      }
      j = next.index + length(as.character(next.val)) + 1
      seed.val <- substring(seed,seed.val.index, seed.val.index)
      seed.string <- paste(seed.string, paste(replicate(next.val,as.character(seed.val)),collapse=""), sep="")
    }
    seed.string
  }
}


# Creates data frame from stats file path, and allows for additional information
# to be added to the data frame.
# Arguments:
#   fname - path to stats file being analyzed
#   include.seeds - if true, read seed file in same dir as fname and put seed info in df
#   include.type  - if true, determine whether seq/sim analysis and include info in df
#   include.chrom - if true, determine what chromosome was analyzed and include in df
#   include.dir   - if true, include base name of directory in which fname is located in df
addInfo <- function(fname, include.seeds, include.type, include.chrom, include.dir){
  dat <- getDf(fname) # read stats file into data frame

  if(include.type) {  # assumes that type information is stored 2 directories above
    type <- basename(dirname(dirname(fname)))
    dat$type <- rep(type,nrow(dat))
  }
  if(include.chrom){  # assumes that chrom information is stored 3 directories above
    chrom <- basename(dirname(dirname(dirname(fname))))
    dat$chrom <- rep(chrom,nrow(dat))
  }
  dir.name <- basename(dirname(fname))
  if(include.dir){
    dat$dir <- rep(dir.name,nrow(dat))
  }

  if(include.seeds){
    seeds <- read.table(file=paste(dirname(fname),"seed_file.txt",sep="/"),colClasses="character")
    mylist <- rep(NA,nrow(dat))
    seed.lengths <- rep(NA,nrow(dat))
    seed.widths <- rep(NA,nrow(dat))
    seed.wlRatios <- rep(NA,nrow(dat))
    seed.strings <- rep(NA,nrow(dat))
    for (i in 1:nrow(dat)){
      s1.index <- dat[i,"seed"]
      s2.index <- match(s1.index,seeds[,1],nomatch=-1)

      if (s2.index != -1){
        seed <- as.character(seeds[s2.index, 2])
        w <- numberOf(seed, 1)
        l <- w + numberOf(seed, 0)

        mylist[i] <- seed
        seed.widths[i] <- w
        seed.lengths[i] <- l
        seed.wlRatios[i] <- w/l
        seed.strings[i] <- toSeedString(seed)
      }
      else{
        mylist[i] <- NA
      }
    }
    names(dat) <- sub("seed", "seed#", names(dat))
    dat$seed <- mylist
    dat$length <- seed.lengths
    dat$width <- seed.widths
    dat$wlRatio <- seed.wlRatios
    dat$seedString <- seed.strings
  }
  freq <- getMatches("(?<=F)[[:digit:]]+",dir.name)
  dat$freq <- rep(freq, nrow(dat))

  if("cc" %in% names(dat)){
    dat$cover <- dat[,"cc"]
  }
  dat
}

# Sorting function. Sorts data from dat according to formula form provided.
# Formula grammar rules:  Use + for ascending, - for descending.
#                         Sorting is left to right in the formula
sortBy <- function(form,dat){
  if(form[[1]] != "~") {
    stop("Formula must be one-sided.")
  }
  # Make the formula into character and remove spaces
  formc <- as.character(form[2]) #substr(form,2,nchar(form)))
  formc <- gsub(" ","",formc)
  # Extract the variables from the formula
  vars <- unlist(strsplit(formc, "[\\+\\-]"))
  vars <- vars[vars!=""] # Remove spurious "" terms
  # Build a list of arguments to pass to "order" function
  calllist <- list()
  pos=1 # Position of + or -
  for(i in 1:length(vars)){
    varsign <- substring(formc,pos,pos)
    pos <- pos+1+nchar(vars[i])
    if(is.factor(dat[,vars[i]])){
      if(varsign=="-")
        calllist[[i]] <- -rank(dat[,vars[i]])
      else
        calllist[[i]] <- rank(dat[,vars[i]])
    }
    else {
      if(varsign=="-")
        calllist[[i]] <- -dat[,vars[i]]
      else
        calllist[[i]] <- dat[,vars[i]]
    }
  }
  dat[do.call("order",calllist),]
}

# Removes "extra" information, as specified. If more_info is not required, removes
# all timing and memory usage information.
removeExcess <-function(dat, more_info){
  drops <- c("cc")
  if(!more_info){
    drops2 <- c("ToolCpuTime","ToolWallTime","ToolMem","ToolVMem","RMCpuTime","RMWallTime","RMMem","RMVMem")
    drops <- c(drops, drops2)
  }
  else{
    #names(dat) <- sub("Tool", "", names(dat))
    #names(dat) <- sub("Cpu", "Cpu\t", names(dat))
    #names(dat) <- sub("RM", "rm", names(dat))
    #names(dat) <- sub("Mem", "Memory", names(dat))
    #names(dat) <- sub("X.tool", "ToolName", names(dat))
    names(dat) <- sub("Time", "", names(dat))
    #names(dat) <- sub("Tool", "", names(dat))
  }
  dat[,!(names(dat) %in% drops)]
}

trim <- function(x) {
    gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

# Formats data from data frame dat, removing any "excess" info and sorting dat
# according to formula provided. Outputs information according to user specification.
# Arguments:
#   dat     - data frame to be formatted/output
#   verbose - if true, formatted data will be output to stdout
#   output  - file to output formatted data to
#   form    - formula specifying sort order of dat
#   auto_out- if true, creates automated name for file and outputs to file
#   dir     - path in which any output files are located
#   more_info- if false, timing/memory usage info not included in output
getSorted <- function(dat, verbose, output, form, auto_out, dir, more_info) {
  dat <- removeExcess(dat, more_info)
  sortedData <- sortBy(formula(form), dat)
  name.width <- nchar("rpt_scout")
  sortedData <-format(sortedData, width = name.width, justify = "left")
  names(sortedData) <- format(names(sortedData), width = name.width, justify = "left")

  if (!is.na(output)) {
      write.table(sortedData, file=output, sep="\t", row.names=FALSE)
  }
  else if (auto_out) {
      file.with.date <- paste("sorted_stats.",trim(as.character(format(Sys.time(), "%Y.%m.%d"))),".txt", sep="")
      write.csv(sortedData, file=paste(dir, file.with.date, sep="/"), row.names=FALSE)
  }
  if (verbose) {
      write.fwf(sortedData, sep=",  ")
  }
}
