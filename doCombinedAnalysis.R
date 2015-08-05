#!/software/R/3.0.2-gcc/bin/Rscript

# doCombinedAnalysis.R : script for combining stats data from multiple RAIDER_eval
#   output directories. Recursively searches for all files in specified directory
#   that are named "stats.txt". Combines all data together, and allows options for
#   adding/removing/formatting data.
#   Arguments:
#     -v/--verbose  specifies whether output will print to stdout. default = True
#     -d/--dir      specifies directory in which to recursively search for stats files
#     -o/--output   specifies file to print output to
#     --formula     specifies sorting formula for output. default = ~-tpr (Descending by sensitivity)
#     -a/--auto_out specifies to print output to default file in same directory
#     -t/--add_type specifies to add info about whether analyzing seq/sim data
#     -c/--add_chrom specifies to add info about chromosome being analyzed
#     -s/--suppress_seeds specifies not to add seed specific information
#     -i            specifies to not include time/space complexity information
#     -m            specifies to include directory information
# by Carly Schaeffer

suppressPackageStartupMessages(require(optparse))
source("sorted_analyze.R")

# allowed arguments
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
}

combineAnalyze <- function(verbose, output, form, auto_out, dir, include.seeds, include.type, include.chrom, more_info, include.dir) {
  files <- list.files(path=dir,pattern="stats.txt",recursive=TRUE,full.names =TRUE)
  files<-files[basename(files) == "stats.txt"]
  myData <- combineData(files, include.seeds, include.type, include.chrom, include.dir)
  source("sorted_analyze.R")
  getSorted(myData, verbose, output, formula(form), auto_out, dir, more_info)
}

combineAnalyze(args2$verbose, args2$output, args2$formula, args2$auto_out, args2$dir, args2$suppress_seeds, args2$add_type, args2$add_chrom, args2$suppress_complexity_info, args2$suppress_dir_info)
