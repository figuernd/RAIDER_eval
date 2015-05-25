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
    make_option(c("--formula"), type="character", default="~+tpr", 
        help = "Sorting formula [default \"%(default)s\"]"),
    make_option(c("-a", "--auto_out"), action="store_true", default=FALSE,
        help="Automatically create output file in same directory as stats file"),
    make_option(c("-t", "--add_type"), action="store_true", default=FALSE,
        help="Automatically add type info (sim vs seq_files)."),
    make_option(c("-c", "--add_chrom"), action="store_true", default=FALSE,
        help="Automatically add chrom info."),
    make_option(c("-s", "--add_seeds"), action="store_true", default=FALSE,
        help="Automatically add seed info."),
    make_option(c("-i", "--more_info"), action="store_true", default=FALSE,
        help="Include info about time and space"),
    make_option(c("-m", "--even_more_info"), action="store_true", default=FALSE,
        help="Include dir name")
)

args = parse_args(OptionParser(option_list=option_list))


sortedAnalyze <- function(fname, verbose, output, form, auto_out, add_seeds, add_type, add_chrom, more_info, even_more_info) {
    if(!is.na(fname)) {
        # getSorted(read.table(fname,header=TRUE,comment.char=""), verbose, output, form, auto_out, dirname(fname))
        df <- addInfo(fname, add_seeds, add_type, add_chrom, even_more_info) #getDf(fname)
        getSorted(df, verbose, output, form, auto_out, dirname(fname), more_info)
    } else {
        cat("Didn't specify file name\n", file=stderr()) # print error messages to stderr
    }
}

sortedAnalyze(args$fname, args$verbose, args$output, args$formula, args$auto_out, args$add_seeds, args$add_type, args$add_chrom, args$more_info, args$even_more_info)
 

