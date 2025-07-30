#!/usr/bin/env Rscript

library(optparse)
library(Rsamtools)
library(data.table)
library(stringr)

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

option_list <- list(
    make_option(c("-i", "--index"), action="store", default=NA, type='character', help="FASTA index file"),
    make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output FASTA file"),
    make_option(c("-f", "--fst-cutoff"), action="store", default=0.7, type='double', help="Fst cutoff value")
)

# use debug arguments if we have 'em
opt_args <- if (exists("debug_args")) debug_args else commandArgs(TRUE)

# parse command-line options
opt <- parse_args(
    OptionParser(
        option_list=option_list,
        formatter=nice_formatter,
        prog="extractsequences.R",
        usage="%prog [options] <fasta_file> <fst_file>"
    ),
    convert_hyphens_to_underscores = TRUE,
    positional_arguments = 2,
    args = opt_args
)

index_file <- opt$options$index
output_file <- opt$options$output
fst_cutoff <- opt$options$fst_cutoff

if (!file.exists(opt$args[1])) {
  stop(sprintf("FASTA file %s does not exist.",opt$args[1]))
}
if (!file.exists(opt$args[2])) {
  stop(sprintf("Fst file %s does not exist.",opt$args[2]))
}
if (!is.na(index_file) & !file.exists(index_file)) {
  stop(sprintf("Index file %s does not exist.",index_file))
}

# load indexed fasta file & get sequence lengths
fasta <- FaFile(file=opt$args[1],index=index_file)
lengths <- seqlengths(fasta)

# use data.table for speed & memory (presumably)
# get sites with "strong" fst and create genomic ranges from them
fst <- fread(opt$args[2])[fst > fst_cutoff, ]
fst <- fst[,.(methods = str_c(unique(method),collapse=',')), by=list(chrom,pos)][order(chrom,pos)]
fst <- fst[,.(positions = str_c(pos,"[",methods,"]",collapse=";")), by = chrom]
fst[, `:=`(len = lengths[chrom],name=str_c(chrom,' positions=',positions))]
fst[,range := str_c(chrom,':1-',len)]

# get pull sequences from fasta using genomic ranges
ranges <- GRanges(fst$range)
seqs <- getSeq(fasta,ranges)
# assign new names to sequences including strongly differentiated positions
names(seqs) <- fst$name

# save to fasta
writeXStringSet(seqs,output_file,format="fasta")
