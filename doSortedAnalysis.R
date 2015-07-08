#!/software/R/3.0.2-gcc/bin/Rscript

suppressPackageStartupMessages(require(optparse))
source("sorted_analyze.R")

# specify our desired options 
# by default ArgumentParser will add an help option 
option_list = list(
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
        help="Print extra output [default]"),
    make_option(c("-f", "--fname"), action="store", type="character", default=NA,
        help="Statistics file name"),
    make_option(c("-o", "--output"), type="character", default=NA, 
        help="File to print output to"),
    make_option(c("--formula"), type="character", default="~-tpr", 
        help = "Sorting formula [default \"%(default)s\"]"),
    make_option(c("-a", "--auto_out"), action="store_true", default=FALSE,
        help="Automatically create output file in same directory as stats file"),
    make_option(c("-t", "--add_type"), action="store_true", default=FALSE,
        help="Automatically add type info (sim vs seq_files)."),
    make_option(c("-c", "--add_chrom"), action="store_true", default=FALSE,
        help="Automatically add chrom info."),
    make_option(c("-s", "--suppress_seeds"), action="store_false", default=TRUE,
        help="Do not include seed info."),
    make_option(c("-i", "--suppress_complexity_info"), action="store_false", default=TRUE,
        help="Do not include info about time and space complexity"),
    make_option(c("-m", "--include_dir_info"), action="store_true", default=FALSE,
        help="Include directory name")
)

args = parse_args(OptionParser(option_list=option_list))


sortedAnalyze <- function(fname, verbose, output, form, auto_out, suppress_seeds, add_type, add_chrom, suppress_complexity_info, include_dir_info) {
    if(!is.na(fname)) {
        # getSorted(read.table(fname,header=TRUE,comment.char=""), verbose, output, form, auto_out, dirname(fname))
        df <- addInfo(fname, suppress_seeds, add_type, add_chrom, include_dir_info) #getDf(fname)
        getSorted(df, verbose, output, form, auto_out, dirname(fname), suppress_complexity_info)
    } else {
        cat("Didn't specify file name\n", file=stderr()) # print error messages to stderr
    }
}

sortedAnalyze(args$fname, args$verbose, args$output, args$formula, args$auto_out, args$suppress_seeds, args$add_type, args$add_chrom, args$suppress_complexity_info, args$include_dir_info)
 

