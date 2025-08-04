#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(poolfstat)
  library(optparse)
  # library(tibble)
  # library(stringr)
  # library(readr)
  # library(dplyr)
  library(data.table)
})


# help message formatter
nice_formatter <- function(object) {
    cat(object@usage, fill = TRUE)
    cat(object@description, fill = TRUE)
    cat("Options:", sep = "\n")

    options_list <- object@options
    for (ii in seq_along(options_list)) {
        option <- options_list[[ii]]
        cat("  ")
        if (!is.na(option@short_flag)) {
            cat(option@short_flag)
            if (optparse:::option_needs_argument(option)) {
                cat(" ", toupper(option@metavar), sep = "")
            }
            cat(", ")
        }
        if (!is.null(option@long_flag)) {
            cat(option@long_flag)
            if (optparse:::option_needs_argument(option)) {
                cat("=", toupper(option@metavar), sep = "")
            }
        }
        cat("\n    ")
        cat(sub("%default", optparse:::as_string(option@default), option@help))
        cat("\n\n")
    }
    cat(object@epilogue, fill = TRUE)
    return(invisible(NULL))
}

list_opt <- function(opt, flag, val, parser, ...) {
    unlist(strsplit(val,","))
}

# set up option list
option_list <- list(
    make_option(c("-w", "--window-size"), action="store", default=0, type='double', help="Fst window size"),
    make_option(c("-c", "--min-coverage"), action="store", default=-1, type='double', help="Minimum coverage per pool"),
    make_option(c("-C", "--min-count"), action="store", default=2, type='double', help="Minimum Minimal allowed read count per base"),
    make_option(c("-x", "--max-coverage"), action="store", default=1e06, type='double', help="Maximum coverage per pool"),
    make_option(c("-m", "--min-maf"), action="store", default=0.01, type='double', help="Minimum minor allele frequnency"),
    make_option(c("-i", "--indels"), action="store_true", default=FALSE, type='logical', help="Keep indel variants"),
    make_option(c("-d", "--headers"), action="store_true", default=FALSE, type='logical', help="Include headers in output files"),
    make_option(c("-p", "--prefix"), action="store", default="", type='character', help="Output file prefix"),
    make_option(c("-s", "--pool-sizes"), action="callback", default=NA, type='character', help="Pool sizes",callback=list_opt),
    make_option(c("-n", "--pool-names"), action="callback", default=NA, type='character', help="Pool names",callback=list_opt),
    make_option(c("-t", "--threads"), action="store", default=1, type='double', help="Number of threads"),
    make_option(c("-l", "--lines-per-thread"), action="store", default=1e06, type='double', help="Number of lines to load per thread")
)

# use debug arguments if we have 'em
if (exists('debug_args')) {
    opt_args <- debug_args
} else {
    opt_args <- commandArgs(TRUE)
}

# parse command-line options
opt <- parse_args(
    OptionParser(
        option_list=option_list,
        formatter=nice_formatter,
        prog="poolfstat.R",
        usage="%prog [options] <sync_file>"
    ),
    convert_hyphens_to_underscores = TRUE,
    positional_arguments = 1,
    args = opt_args
)


# load pool data from sync
pooldata <- popsync2pooldata(
    sync.file = opt$args[1],
    poolsizes = as.numeric(opt$options$pool_sizes),
    poolnames = opt$options$pool_names,
    min.rc = opt$options$min_count,
    min.cov.per.pool = opt$options$min_coverage,
    max.cov.per.pool = opt$options$max_coverage,
    min.maf = opt$options$min_maf,
    noindel = !opt$options$indels,
    nthreads = opt$options$threads
)

# calculate F-stats
fst <- computeFST(pooldata,sliding.window.size=opt$options$window_size)

# get population combo string
# fn <- str_c(pooldata@poolnames,collapse='-')
fn <- paste0(pooldata@poolnames,collapse='-')

# save global Fst
# global_fst <- tibble(pop1=pooldata@poolnames[1],pop2=pooldata@poolnames[2],fst=fst$Fst[1])
# write_tsv(global_fst,str_glue("{opt$options$prefix}_global_{fn}.tsv"))
global_fst <- data.table(pop1=pooldata@poolnames[1],pop2=pooldata@poolnames[2],fst=fst$Fst[1])
fwrite(global_fst,sprintf("%s_global_%s.tsv",opt$options$prefix,fn),sep="\t")

# save sliding window Fst, if they exist
if (!is.null(fst$sliding.windows.fvalues)) {
    # sliding <- as_tibble(fst$sliding.windows.fvalues) %>%
    #     rename(chrom=1,start=2,end=3,mid=4,cum_mid=5,fst=6)
    # write_tsv(sliding,str_glue("{opt$options$prefix}_sliding-{opt$options$window_size}_{fn}.tsv"))
  setDT(fst$sliding.windows.fvalues)
  setnames(fst$sliding.windows.fvalues,new=c('chrom','start','end','mid','cum_mid','fst'))
  fwrite(fst$sliding.windows.fvalues,sprintf("%s_sliding-%d_%s.tsv",opt$options$prefix,opt$options$window_size,fn),sep="\t")
}

# pools <- str_c(str_c(pooldata@poolnames,collapse=':'),".fst")
pools <- paste0(paste0(pooldata@poolnames,collapse=':'),".fst")

# save per-snp Fst, match popoolation output
# snp_fst <- pooldata@snp.info %>%
#     as_tibble() %>%
#     select(chrom=1,pos=2) %>%
#     mutate(
#       window_size = opt$options$window_size,
#       covered_fraction = 1,
#       # depth = rowSums(pooldata@readcoverage),
#       avg_min_cov = do.call(pmin,as.data.frame(pooldata@readcoverage)),
#       pop1 = opt$options$pool_names[1],
#       pop2 = opt$options$pool_names[2],
#       fst = fst$snp.Fstats$Fst
#     ) %>%
#     select(chrom,pos,window_size,covered_fraction,avg_min_cov,pop1,pop2,fst)
# write_tsv(snp_fst,str_glue("{opt$options$prefix}_{fn}.fst"),col_names=opt$options$headers)

setDT(pooldata@snp.info)
cols <- names(pooldata@snp.info)
pooldata@snp.info[,cols[-c(1:2)] := NULL]
setnames(pooldata@snp.info,c('chrom','pos'))
pooldata@snp.info[,`:=`(
  window_size = opt$options$window_size,
  covered_fraction = 1,
  avg_min_cov = do.call(pmin,as.data.table(pooldata@readcoverage)),
  pop1 = opt$options$pool_names[1],
  pop2 = opt$options$pool_names[2],
  fst = fst$snp.Fstats$Fst
)]
fwrite(pooldata@snp.info,sprintf("%s_%s.fst",opt$options$prefix,fn),col.names = opt$options$headers)
