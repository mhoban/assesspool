#!/usr/bin/env Rscript

# load packages ---------------------------------------------------------------------------------------------------

suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(fs))


# setup -----------------------------------------------------------------------------------------------------------
"
Usage: prep_sync.R [options] <variant> <genome_ref>

Options:
  -p --prefix <prefix> Output file prefix [default: sync]
  -m --multiallelic    Keep multiallelic SNPs
  -d --indels          Keep indel loci
  -s --split           Split resulting sync file by population pairs
  -r --header          Include header in sync file (for grenedalf)
  -c --combined        Save combined sync file
  -h --help            This message
" -> doc

if (exists('debug_args')) {
    args <- debug_args
} else {
    args <- commandArgs(TRUE)
}

# process command-line options
opts <- docopt(doc,args)

if (!file_exists(opts$variant)) {
    stop(str_glue("VCF file {opts$variant} does not exist!"))
}

if (!file_exists(opts$genome_ref)) {
    stop(str_glue("Reference file {opts$genome_ref} does not exist!"))
}

if (!opts$split & !opts$combined) {
    stop(str_glue("One or more of --split or --combined must be included!"))
}

# read a vcf file with the specified fields
read_vcf <- function(
        vcf_file, ref_file, fixed = c("ALT","QUAL"),
        info = c("NS","DP","TYPE","LEN","NUMALT"),
        per_sample = c("DP","RO","AO","GT","GQ")) {
    param <- ScanVcfParam(fixed = fixed, info = info, geno = per_sample)
    suppressWarnings( readVcf(vcf_file, ref_file, param) )
}

# add column(s) to a tibble if they're missing
# (swiped and modified from poissonconsulting/tidyplus)
add_missing_column <- function(
        .data, ..., .before = NULL, .after = NULL,
        .name_repair = c("check_unique", "unique", "universal", "minimal")) {
    if (!is.data.frame(.data)) {
        stop("add_column(.data = 'must be a data frame')")
    }
    dots <- rlang::list2(...)
    if (length(dots) == 0L) {
        return(.data)
    }
    names <- names(dots)
    names <- names[!names %in% colnames(.data)]
    if (!length(names)) {
        return(.data)
    }
    tibble::add_column(.data, !!!dots[names],
                       .before = .before, .after = .after,
                       .name_repair = .name_repair
    )
}


# load vcf and populate tibbles -----------------------------------------------------------------------------------

# columns of interest:
# DP: Total read depth at the locus
# RO: Count of full observations of the reference haplotype
# AO: Count of full observations of the alternate haplotype
vcf_cols <- c("DP","RO","AO")

# load vcf from file
vcf <- read_vcf(opts$variant,opts$genome_ref)

# get sample/population names
POPS <- samples(header(vcf))

# create 'header' tibble
vcf_info <- info(vcf) %>%
    # cast as tibble with rownames as snp ID
    as_tibble(rownames="snpid") %>%
    clean_names() %>%
    # break out chromosome and position from snp ID
    separate_wider_regex(snpid,pattern=c(chrom='.+',':',pos='.+','_.*$'),cols_remove = FALSE) %>%
    # put snp ID first
    dplyr::select(snpid,everything()) %>%
    dplyr::rename(tdp=dp,an=numalt) %>%
    mutate(
        # get number of alternate alleles (I think)
        alternates = map_int(type,length),
        # make pos column numeric,
        pos = as.numeric(pos),
        # get reference allele
        ref = as.character(rowRanges(vcf)$REF),
        # get list of alternate alleles
        alt = as.list(CharacterList(rowRanges(vcf)$ALT)),
        # length of alternate allele(s)
        alen = map(alt,str_length),
        # length of reference allele
        rlen = str_length(ref),
        # length of insert (if any)
        ins_len = map2(alen,rlen,~.x-.y),
        # length of deletion (if any)
        del_len = map2(rlen,alen,~.x-.y),
        # set type to first type in each list
        type = map_chr(type,dplyr::first)
    )


# get all sample/header data in long format
vcf_data_long <- geno(vcf)[vcf_cols] %>%
    # map through depth variables
    map(\(data) {
        # get counts
        data %>%
            # cast as tibble with rownames as snp ID
            as_tibble(rownames="snpid") %>%
            mutate(
                # replace NAs with zeroes
                across(all_of(POPS),\(x) {
                    if (is.list(x)) map(x,\(y) replace_na(y,0))
                    else replace_na(x,0)
                })
            ) %>%
            # get rid of 'group' columns
            dplyr::select(-starts_with("group"))
    }) %>%
    # map through datasets
    imap(\(data,var_name) {
        data %>%
            # pivot to longer version
            pivot_longer(all_of(POPS),names_to="pop",values_to = var_name)
    }) %>%
    # smash them all together with full joins
    purrr::reduce(\(a,b) full_join(a,b,by=c("snpid","pop"))) %>%
    # order by chrom, pos, and population
    # arrange(chrom,pos,pop) %>%
    # join vcf info dataset
    full_join(vcf_info,by="snpid") %>%
    clean_names() %>%
    # order columns
    dplyr::select(snpid,chrom,pos,alternates,type,an,rlen,alen,ins_len,del_len,ref,ro,alt,ao,pop) %>%
    # sort by snp ID, chromosome, position
    arrange(snpid,chrom,pos)

# clear up memory
gc()

# generate sync file(s) -------------------------------------------------------------------------------------------

# sync file columns
bases <- c("A","T","C","G","N","DEL")
# named list of zeroes to use with add_missing_column
base_cols <- bases %>%
    set_names() %>%
    map(~0)

# generate sync table
sync <- vcf_data_long %>%
    dplyr::filter(
        # include_mulltiallelic | (alternates == 1 & an == 1 & rlen == 1),
        opts$multiallelic | alternates == 1,
        opts$indels | !(type %in% c('ins','del')),
        type != "snp" | (ref %in% bases | alt %in% bases)
    )    %>%
    # unnest list columns (expands multiallelic snps)
    unnest(where(is.list)) %>%
    # make sure we're only dealing with single-base insertions/deletions
    dplyr::filter(!(type %in% c('ins','del')) | (ins_len == 1 | del_len == 1)) %>%
    # get ref/alt counts and bases
    mutate(
        # get reference base
        refbase = case_when(
            # for deletions, the second to last base of the reference
            type == "del" ~ str_sub(ref,start=rlen-del_len,end=rlen-1),
            # for deletions, the second to last base of the alternate
            type == "ins" ~ str_sub(alt,start=alen-ins_len,end=alen-1),
            # otherwise just use ref
            .default = ref
        ),
        # rename base as DEL for indels
        altbase = case_when(
            type %in% c('ins','del') ~ "DEL",
            .default = alt
        ),
        # get correct reference count
        refcount = case_when(
            type == "ins" ~ ao,
            .default = ro
        ),
        # get correct alternate count
        altcount = case_when(
            type == "ins" ~ ro,
            .default = ao
        )
    ) %>%
    # save reference base
    mutate(ref = refbase) %>%
    # filter to single-base snps
    filter(refbase %in% bases & altbase %in% bases) %>%
    # keep columns
    dplyr::select(snpid,chrom,pos,pop,alternates,type,refbase,altbase,refcount,altcount,ref) %>%
    # pivot ref bases/counts to wide
    pivot_wider(names_from = "refbase", values_from = "refcount", values_fill = 0) %>%
    # pivot alt bases/counts to wide, duplicate names will look like A--->2, G--->2, etc.
    pivot_wider(names_from = "altbase", values_from = "altcount", values_fill = 0, names_repair = \(names) make.unique(names,sep = "--->") ) %>%
    # pivot base columns back to longer to combine them
    pivot_longer(-c(snpid:ref),names_to = "base", values_to = "count") %>%
    # fix "repaired" base names
    mutate(base = str_replace(base,"--->.+$","")) %>%
    # get counts per SNP (combines multiallelic loci)
    group_by(across(-count)) %>%
    summarise(count = sum(count)) %>%
    ungroup() %>%
    # go back to wide one more time
    pivot_wider(names_from = base, values_from = count) %>%
    add_missing_column(!!!base_cols) %>%
    mutate(variants = alternates+1) %>%
    dplyr::select(snpid,chrom,pos,pop,type,variants,ref,all_of(bases)) %>%
    arrange(snpid,pop) %>%
    unite("sync",all_of(bases),sep=":") %>%
    pivot_wider(names_from="pop",values_from="sync") %>%
    arrange(chrom,pos) %>%
    rename(`#chr`=chrom) %>%
    select(`#chr`,pos,ref,all_of(POPS))

# clear up memory
gc()

# save output files -----------------------------------------------------------------------------------------------


# save all-population sync file
if (opts$combined) {
    dir_create("combined")
    write_tsv(sync,str_glue("combined/{opts$prefix}_all_pops.sync"),col_names = opts$header)
}

# if desired, split sync files
if (opts$split) {
    dir_create("split")

    # get population combinations
    # and save to individual sync files
    POPS %>%
        combn(2) %>%
        array_branch(2) %>%
        walk(\(combo) {
            combo <- sort(combo)
            fn <- str_glue("split/{opts$prefix}_{combo[1]}-{combo[2]}.sync")
            # filter sync by combination and save it
            sync %>%
                select(`#chr`,pos,ref,all_of(combo)) %>%
                write_tsv(fn,col_names = FALSE)
        })
}

