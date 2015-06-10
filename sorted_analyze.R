#!/software/R/3.0.2-gcc/bin/Rscript

suppressPackageStartupMessages(require(optparse))

suppressPackageStartupMessages(require(gdata))
getDf <- function(fname) {
    #dat <- readLines(fname,n=-1)
    #names <- unlist(strsplit(dat[1]," "))
    #names <- names[names != ""]
    #dat <- dat[-1]
    #dat <- dat[grep("INCOMPLETE",dat,invert=TRUE)]
    #dat <- unlist(lapply(dat,strsplit," "))
    #dat <- dat[dat != ""]
    #dat <- data.frame(matrix(dat,ncol=length(names),byrow=TRUE))
    #colnames(dat)<-names
    #as.data.frame(sapply(dat, function(x) gsub("\"", "", x)))
    #dat[2:length(names)] <- lapply(dat[2:length(names)], function(x) as.numeric(levels(x)[x]))
    #dat
    dat <- read.table(fname,header=TRUE,comment.char="",fill=TRUE)
    dat <- dat[rowSums(is.na(dat)) < 4, ]

    #cat(dim(dat))
    #cat(length(names(dat)))
    #write.table(dat)
}


sumOverMatches <- function(matches){
    sum(as.numeric(unlist(matches)))
}

getMatches <- function(expr, seed){
    regmatches(seed, gregexpr(expr, seed, perl=TRUE))
}

numberOf <- function(seed, val){
    match.string <- paste("(?<=", val, "\\^\\{)[[:digit:]]+(?=\\})", sep="")
    numOnes <- sumOverMatches(getMatches(match.string, seed))
    match.string <- paste(val, "(?!\\^)", sep="")
    numOnes + length(unlist(getMatches(match.string, seed)))
}

toSeedString <- function(seed){
    seed.string = ""
    match.string <- "(?<=\\^\\{)[[:digit:]]+(?=\\})"
    match.indices <- gregexpr(match.string,seed,perl=TRUE)[[1]]
    match.values <- regmatches(seed,gregexpr(match.string,seed,perl=TRUE))[[1]]
    j = 1
    for( i in 1:length(match.indices)){
        next.index <- match.indices[i]
        next.val <- match.values[i]
        seed.val.index = next.index - 3
        if(seed.val.index > j + 1){
            print( substring(seed, j, seed.val.index - 1))
            seed.string <- paste(seed.string, substring(seed, j, seed.val.index - 1), sep="")
        }
        j = next.index + length(as.character(next.val)) + 1
        seed.val <- substring(seed,seed.val.index, seed.val.index)
        seed.string <- paste(seed.string, paste(replicate(next.val,as.character(seed.val)),collapse=""), sep="")
    }
    seed.string
    

}


addInfo <- function(fname, include.seeds, include.type, include.chrom, include.dir){
    dat <- getDf(fname) #read.table(fname,header=TRUE,comment.char="")
    if(include.type) {
        type <- basename(dirname(dirname(dirname(fname))))
        dat$type <- rep(type,nrow(dat))
    }
    if(include.chrom){
        chrom <- basename(dirname(dirname(fname)))
        dat$chrom <- rep(chrom,nrow(dat))
    }
    dir.name <- basename(dirname(fname))
    
    if(include.dir){
        dat$dir <- rep(dir.name,nrow(dat))
    }
    
    
    if(include.seeds){
        seeds <- read.table(file=paste(dirname(fname),"seed_file.txt",sep="/"))
        mylist <- rep(NA,nrow(dat)) #list()
        seed.lengths <- rep(NA,nrow(dat))
        seed.widths <- rep(NA,nrow(dat))
        seed.wlRatios <- rep(NA,nrow(dat))
        seed.strings <- rep(NA,nrow(dat))
        for (i in 1:nrow(dat)){
            s1.index <- dat[i,"seed"]
            #print(s1.index)
            #print("\n")
            s2.index <- match(s1.index,seeds[,1],nomatch=-1)
            #print(s2.index)
            #print("\n")
            if (s2.index != -1){
                #print(seeds[s2.index, 2])
                #print("\n")
                seed <- as.character(seeds[s2.index, 2]) #as.character(seeds[s2.index, 2])
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
        #drops <- c("seed")
        #dat <- dat[,!(names(dat) %in% drops)]
    }
    freq <- getMatches("(?<=F)[[:digit:]]+",dir.name)
    dat$freq <- rep(freq, nrow(dat))
    if("cc" %in% names(dat)){
        dat$cover <- dat[,"cc"]
    }
    dat
}

sortBy <- function(form,dat){
    # Use + for ascending, - for decending.
    # Sorting is left to right in the formula
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



removeExcess <-function(dat, more_info){
    drops <- c("cc")#c("tp","tn","fp","fn")
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

getSorted <- function(dat, verbose, output, form, auto_out, dir, more_info) {
    dat <- removeExcess(dat, more_info)
    #return dat
    sortedData <- sortBy(formula(form), dat)
    name.width <- nchar("rpt_scout")# max(sapply(names(sortedData), nchar))
    sortedData <-format(sortedData, width = name.width, justify = "left")
    #print(sortedData[sapply(sortedData, function(x) is.numeric(x))])
    names(sortedData) <- format(names(sortedData), width = name.width, justify = "left")
    #sortedData[sapply(sortedData, function(x) is.numeric(x))] <- lapply(sortedData[sapply(sortedData, function(x) is.numeric(x))], function(x) prettyNum(x))#formatC(x, format="fg", digits=4, flag="#"))
    #print(dat)
    if (!is.na(output)) {
        write.table(sortedData, file=output, sep="\t", row.names=FALSE)
    }
    else if (auto_out) {
        file.with.date <- paste("sorted_stats.",trim(as.character(format(Sys.time(), "%Y.%m.%d"))),".txt", sep="")
        #write.table(sortedData, file=paste(dir,"sorted_stats.txt",sep="/"), sep="\t", row.names=FALSE)
        #write.table(sortedData, file=paste(dir, file.with.date, sep="/"), sep="\t", row.names=FALSE)
    
        #if (and_csv) {
        write.csv(sortedData, file=paste(dir, file.with.date, sep="/"), row.names=FALSE)
        #}
    }
    if (verbose) {
        #write.fwf(sortedData, sep="  ")
        #write.table(sortedData, sep="\t", row.names=FALSE)
        write.fwf(sortedData, sep=",  ")
    }
    #cat(capture.output(sortedData), file = 'dframe.txt', sep = '\n')
    #print(sortedData)
}

 

