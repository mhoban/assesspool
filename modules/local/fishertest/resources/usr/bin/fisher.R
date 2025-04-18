#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(matrixStats))

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

# geometric mean
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}


pval_callback <- function(opt, flag, val, parser, ...) {
  val <- char.expand(val,c("multiply","geometric-mean"))
  switch(
    val,
    'multiply' = prod,
    'geometric-mean' = gm_mean
  )
}

expand_callback <- function(opt, flag, val, parser, args, ...) {
  char.expand(val,args)
}

pops_callback <- function(opt, flag, val, parser, ...) {
  if (file.exists(val)) {
    read_lines(file)
  } else {
    unlist(strsplit(val,","))
  }
}


# set up option list
option_list <- list(
    make_option(c("-p", "--pools"), action="callback", default=NA, type='character',
                help="Subset of pools for which to calculate Fisher test (comma-separated or listed in a text file)",callback=pops_callback),
    make_option(c("-m", "--window-type"), action="callback", default="single", type='character', help="Window type",
                callback=expand_callback, callback_args = list(args=c("single","interval","queue"))),
    make_option(c("-w", "--window-size"), action="store", default=1, type='double', help="Fst window size"),
    make_option(c("-s", "--window-step"), action="store", default=0, type='double', help="Fst window step"),
    make_option(c("-C", "--min-count"), action="store", default=2, type='double', help="Minimum Minimal allowed read count per base"),
    make_option(c("-c", "--min-coverage"), action="store", default=0, type='double', help="Minimum coverage per pool"),
    make_option(c("-x", "--max-coverage"), action="store", default=1e03, type='double', help="Maximum coverage per pool"),
    make_option(c("-a", "--pval-aggregate"), action="callback", default=prod, type='character', help="P-value aggregation method",callback=pval_callback),
    make_option(c("-A", "--all-snps"), action="store_true", default=FALSE, type='logical', help="Save fisher test results for all SNPs, regardless of window size"),
    make_option(c("-j", "--adjust-pval"), action="store_true", default=FALSE, type='logical', help="Adjust p-value for multiple comparisons"),
    make_option(c("-J", "--adjust-method"), action="store", default="bonferroni", type='character', help="P-value adjustment method"),
    make_option(c("-f", "--prefix"), action="store", default="", type='character', help="Output file prefix"),
    make_option(c("-o", "--output-dir"), action="store", default=".", type='character', help="Output directory"),
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
        prog="fisher.R",
        usage="%prog [options] <frequency_file>"
    ),
    convert_hyphens_to_underscores = TRUE,
    positional_arguments = 1,
    args = opt_args
)

window_type  <- opt$options$window_type
if (window_type == "single") {
  opt$options$window_size <- 1
  opt$options$window_step <- 1
}

window_size   <- opt$options$window_size
window_step   <- (if (opt$options$window_step) opt$options$window_step else window_size)
min_count     <- opt$options$min_count
min_depth     <- opt$options$min_coverage
threads       <- opt$options$threads
thread_lines  <- opt$options$lines_per_thread
file_prefix   <- opt$options$prefix
outdir        <- opt$options$output_dir
p_combine     <- opt$options$pval_aggregate
save_all      <- opt$options$all_snps
adjust_p      <- opt$options$adjust_pval
adjust_method <- opt$options$adjust_method
pools          <- opt$options$pools
if (!all(is.na(pools))) {
  if (length(pools) != 2) {
    stop("--pools must be a list of exactly two pool names, either comma-separated or in a text file")
  }
}

if (!dir.exists(outdir)) {
  stop("Output directory does not exist.")
}

if (!file.exists(opt$args[1])) {
  stop(sprintf("Frequency file %s does not exist.",opt$args[1]))
}

# filter frequency table down to just what we'd consider snps
# this filtering should match popoolation/poolfst

# split(ceiling(seq_along(.)/50)) %>%
  
freq_snps <- read_csv(opt$args[1],col_types = cols(),progress = FALSE) 

if (all(is.na(pools))) {
  pools <- freq_snps %>%
    slice(1) %>%
    pivot_longer(-c(CHROM:ALT,starts_with("TOTAL",ignore.case = FALSE)),names_to = c("pool","measure"),values_to="count",names_pattern = "^(.+)\\.([^.]+)$") %>%
    pull(pool) %>%
    sort() %>%
    unique()
}

# make sure we're dealing with two pools
npool <- length(pools)
# and bail if we're not
if (npool != 2) {
  stop("Frequency file must contain exactly two pools")
}

freq_snp_og <- freq_snps
# continue filtering
freq_snps <- freq_snp_og %>%
  select(CHROM:ALT,starts_with(paste0(pools,"."),ignore.case = FALSE)) %>%
  mutate(
    `TOTAL.REF_CNT` = rowSums(pick(ends_with(".REF_CNT",ignore.case = FALSE))),
    `TOTAL.ALT_CNT` = rowSums(pick(ends_with(".ALT_CNT",ignore.case = FALSE))),
    `TOTAL.DEPTH` = rowSums(pick(ends_with(".DEPTH",ignore.case = FALSE))),
  ) %>%
  mutate(
    lwr = floor((POS-1)/window_step)*window_step+1,
    upr = lwr+window_size-1,
    middle=floor((upr+lwr)/2)
  ) %>%
  filter( if_all(ends_with(".DEPTH"),~.x >= min_depth) ) %>%
  add_count(CHROM,middle,name="snp_count") %>%
  rename_with(make_clean_names,.cols = starts_with("TOTAL.",ignore.case = FALSE)) %>%
  filter(
    if_all(starts_with("total_",ignore.case = FALSE) & ends_with("_cnt",ignore.case = FALSE),~.x >= min_count),
    total_ref_cnt != total_depth,
    total_depth >= min_depth
  ) %>%
  rename_with(CHROM:ALT,.fn=make_clean_names) %>%
  select(chrom,pos,middle,ref:alt,ends_with(".REF_CNT"),ends_with(".ALT_CNT"),starts_with("total_"),everything())

fisher_results <- freq_snps %>%
  select(order(colnames(.))) %>%
  select(
    chrom, pos, middle, snp_count,start=lwr,end=upr,
    ends_with(".DEPTH",ignore.case = FALSE),
    ends_with(".REF_CNT",ignore.case = FALSE),
    ends_with(".ALT_CNT",ignore.case = FALSE)
  ) %>%
  pivot_longer(
    -c(chrom,pos,middle,snp_count,start,end,ends_with(".DEPTH",ignore.case = FALSE)),
    names_to = c("pop","measure"),
    values_to = "count",
    names_pattern = "^(.+)\\.([^.]+)$"
  ) %>%
  group_by(chrom,pos,middle,start,end,snp_count,across(ends_with(".DEPTH",ignore.case = FALSE))) %>%
  summarise(pval = fisher.test(matrix(count,ncol=2))$p.value) %>%
  ungroup() %>%
  mutate(min_cov = matrixStats::rowMins(as.matrix(pick(ends_with(".DEPTH",ignore.case = FALSE)))))

if (adjust_p) {
  fisher_results <- fisher_results %>%
    mutate(p_adj = p.adjust(pval,method=adjust_method))
}


pool_col <- paste0(paste0(pools,collapse=":"),".fisher")

if (save_all & window_size > 1) {
  ss <- sprintf("%s/%s_%s_all_snps_fisher.tsv",outdir,file_prefix,paste0(pools,collapse="-"))
  fisher_results %>%
    mutate(
      "{pool_col}" := -log10(pval),
      snps = 1,
      covered_fraction = 1,
      start=pos,
      end=pos,
    ) %>%
    select(chrom,start,middle=pos,end,snps,covered_fraction,avg_min_coverage = min_cov,all_of(pool_col)) %>%
    write_tsv(ss)
}

fisher_results <- fisher_results %>%
  group_by(chrom,middle,start,end,snp_count) %>%
  summarise(
    "{pool_col}" := -log10(p_combine(pval)),
    snps = n(),
    covered_fraction = unique(snp_count)/window_size,
    avg_min_coverage = mean(min_cov)
  ) %>%
  select(chrom,start,middle,end,snps,covered_fraction,avg_min_coverage,all_of(pool_col))  %>%
  ungroup()

ss <- sprintf("%s/%s_%s_window_%d_%d_fisher.tsv",outdir,file_prefix,paste0(pools,collapse="-"),window_size,window_step)
write_tsv(fisher_results,ss)
# hello