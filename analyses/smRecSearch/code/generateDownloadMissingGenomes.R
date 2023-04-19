#!/usr/bin/env Rscript

library(tidyverse)

combine <- function(..., prefix = "", sep = "_") {
  paste0(prefix, levels(interaction(..., sep = sep)))
}

args <- commandArgs(trailingOnly = TRUE)
genomes <- read_lines(args[1]) %>% set_names(.,.)

twobits <- genomes %>% paste0(., ".2bit") %>% set_names(., genomes)
chromsize <- genomes %>% paste0(., ".chrom.sizes") %>% set_names(., genomes)
fasta <- genomes %>% paste0(., ".fa") %>% set_names(., genomes)

twobits %>%
  paste0("rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/", genomes, "/bigZips/", ., " data/2bit/") %>%
  c("#!/bin/bash\n\n", "#Generated using 'generateDownloadMissingGenomes.R'\n\n", .) %>%
  write_lines("code/downloadMissing2bits.bash")

chromsize %>%
  paste0("rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/", genomes, "/bigZips/", ., " data/genome/") %>%
  c("#!/bin/bash\n\n", "#Generated using 'generateDownloadMissingGenomes.R'\n\n", .) %>%
  write_lines("code/downloadMissingchromsize.bash")

fasta %>%
  paste0("rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/", genomes, "/bigZips/", ., " data/genome/") %>%
  c("#!/bin/bash\n\n", "#Generated using 'generateDownloadMissingGenomes.R'\n\n", .) %>%
  write_lines("code/downloadMissingGenome.bash")
