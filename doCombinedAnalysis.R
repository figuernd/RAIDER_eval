#!/software/R/3.0.2-gcc/bin/Rscript
suppressPackageStartupMessages(require(optparse))
source("sorted_analyze.R")

# specify our desired options 
# by default ArgumentParser will add an help option 
option_list_2 = list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
        help="Print extra output [default]"),
    make_option(c("-d", "--dir"), action="store", type="character", default=NA,
        help="Directory to recursively search for stats files"),
    make_option(c("-o", "--output"), type="character", default=NA, 
        help="File to print output to"),
    make_option(c("--formula"), type="character", default="~-tpr", 
        help = "Sorting formula [default \"%(default)s\"]"),
    make_option(c("-a", "--auto_out"), action="store_true", default=FALSE,
        help="Automatically create output file in same directory."),
    make_option(c("-t", "--add_type"), action="store_true", default=FALSE,
        help="Automatically add type info (sim vs seq_files)."),
    make_option(c("-c", "--add_chrom"), action="store_true", default=FALSE,
        help="Automatically add chrom info."),
    make_option(c("-s", "--suppress_seeds"), action="store_false", default=TRUE,
        help="Do not include seed info."),
    make_option(c("-i", "--suppress_complexity_info"), action="store_false", default=TRUE,
        help="Do not include info about time and space complexity"),
    make_option(c("-m", "--suppress_dir_info"), action="store_true", default=FALSE,
        help="Include directory name")
)

args2 = parse_args(OptionParser(option_list=option_list_2))




combineData <- function(files, include.seeds, include.type, include.chrom, include.dir){
    dfs <- lapply(files, function(x) addInfo(x, include.seeds, include.type, include.chrom, include.dir))
    merged <- Reduce(function(x, y) rbind(x,y), dfs)#merge(x, y, all=TRUE), dfs) 
    

    #merged <- merged[complete.cases(merged), ]
}

combineAnalyze <- function(verbose, output, form, auto_out, dir, include.seeds, include.type, include.chrom, more_info, include.dir) {
    files <- list.files(path=dir,pattern="stats.txt",recursive=TRUE,full.names =TRUE)
    files<-files[basename(files) == "stats.txt"]
    myData <- combineData(files, include.seeds, include.type, include.chrom, include.dir)
    source("sorted_analyze.R")
    getSorted(myData, verbose, output, formula(form), auto_out, dir, more_info)
}

combineAnalyze(args2$verbose, args2$output, args2$formula, args2$auto_out, args2$dir, args2$suppress_seeds, args2$add_type, args2$add_chrom, args2$suppress_complexity_info, args2$suppress_dir_info)
 

