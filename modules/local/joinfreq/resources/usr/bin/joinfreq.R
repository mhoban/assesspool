#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr) 
  library(tidyr) 
  library(readr)
  library(janitor) 
  library(purrr)
  library(optparse) 
  library(stringr) 
})

options(dplyr.summarise.inform = FALSE)


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

pops_callback <- function(opt, flag, val, parser, ...) {
  if (file.exists(val)) {
    read_lines(file)
  } else {
    unlist(strsplit(val,","))
  }
}

expand_callback <- function(opt, flag, val, parser, args, ...) {
  char.expand(val,args)
}

option_list <- list(
    make_option(c("-m", "--window-type"), action="callback", default="single", type='character', help="Window type",
                callback=expand_callback, callback_args = list(args=c("single","interval","queue"))),
    make_option(c("-w", "--window-size"), action="store", default=1, type='double', help="Fst window size"),
    make_option(c("-s", "--window-step"), action="store", default=0, type='double', help="Fst window step"),
    make_option(c("-C", "--min-count"), action="store", default=2, type='double', help="Minimum Minimal allowed read count per base"),
    make_option(c("-c", "--min-coverage"), action="store", default=0, type='double', help="Minimum coverage per pool"),
    make_option(c("-x", "--max-coverage"), action="store", default=1e03, type='double', help="Maximum coverage per pool"),
    make_option(c("-f", "--prefix"), action="store", default="", type='character', help="Output file prefix"),
    make_option(c("-o", "--output-dir"), action="store", default=".", type='character', help="Output directory")
)

debug_args <- c(
  "--window-size", "1", "--window-step", "1", "--min-count",
  "2", "--min-coverage", "10", "--max-coverage", "1000", "--prefix",
  "fun", "/Users/deadbilly/projects/pipeline/pool/nf-core-assesspool/testing/test/datatable/testz0r_frequency.csv",
  "/Users/deadbilly/projects/pipeline/pool/nf-core-assesspool/testing/test/datatable/testz0r_fst.csv"
)

# use debug arguments if we have 'em
opt_args <- if (exists("debug_args")) debug_args else commandArgs(TRUE)


# parse command-line options
opt <- parse_args(
    OptionParser(
        option_list=option_list,
        formatter=nice_formatter,
        prog="fisher.R",
        usage="%prog [options] <frequency_file> <fst_file>"
    ),
    convert_hyphens_to_underscores = TRUE,
    positional_arguments = 2,
    args = opt_args
)

window_size   <- opt$options$window_size
window_step   <- (if (opt$options$window_step) opt$options$window_step else window_size)
min_count     <- opt$options$min_count
min_depth     <- opt$options$min_coverage
file_prefix   <- opt$options$prefix
outdir        <- opt$options$output_dir

if (!dir.exists(outdir)) {
  stop("Output directory does not exist.")
}

if (!file.exists(opt$args[1])) {
  stop(sprintf("Frequency file %s does not exist.",opt$args[1]))
}

if (!file.exists(opt$args[2])) {
  stop(sprintf("Fst file %s does not exist.",opt$args[1]))
}

freq_snps <- read_csv(opt$args[1],col_types = cols(),progress = FALSE) 
fst_wide  <- read_csv(opt$args[2],col_types = cols(),progress = FALSE) 

# load fst file, pivot long, and filter out NAs
fst <- fst_wide %>%
  select(chrom:end,ends_with(".fst",ignore.case = FALSE)) %>%
  pivot_longer(-c(chrom:end),names_to=c("pop1","pop2",".value"),names_pattern = "^([^:]+):(.+)\\.(fst)$") %>%
  filter(!is.na(fst),start == end)

freq_snps <- freq_snps %>%
  pivot_longer(-c(CHROM:ALT,starts_with("TOTAL",ignore.case = FALSE)),names_to = c("pool","measure"),values_to="count",names_pattern = "^(.+)\\.([^.]+)$") %>%
  select(CHROM,POS,pool,measure,count) %>%
  janitor::clean_names() %>%
  pivot_wider(names_from="measure",values_from = "count")

# continue filtering
combined <- fst %>%
  left_join(freq_snps,by=c("chrom","start" = "pos","pop1" = "pool")) %>%
  left_join(freq_snps,by=c("chrom","start" = "pos","pop2" = "pool"),suffix=c(".pop1",".pop2")) %>%
  mutate(
    # summarize read depths
    TOTAL.REF_CNT = REF_CNT.pop1 + REF_CNT.pop2,
    TOTAL.ALT_CNT = ALT_CNT.pop1 + ALT_CNT.pop2,
    TOTAL.DEPTH = DEPTH.pop1 + DEPTH.pop2
  ) %>%
  filter(
    # filter by minimum minimum depth and others
    DEPTH.pop1 >= min_depth & DEPTH.pop2 >= min_depth,
    TOTAL.ALT_CNT >= min_count & TOTAL.REF_CNT >= min_count,
    TOTAL.REF_CNT != TOTAL.DEPTH,
    TOTAL.DEPTH >= min_depth
  ) %>%
  mutate( avg_min_cov = pmin(DEPTH.pop1, DEPTH.pop2) ) %>%
  janitor::clean_names() %>%
  arrange(chrom,start,pop1,pop2) 
  
outfile <- sprintf("%s/%s_fst_frequency.tsv",outdir,file_prefix)
write_tsv(combined,outfile)