#!/usr/bin/env Rscript

library(tidyverse)

pt <- read_csv("./data/portTable.csv", comment = '"#"', col_names = c("Translated_BLAT_Port", "Untranslated_BLAT_Port", "Genome", "Name", "Species")) %>%
  filter(!str_detect(Genome, "softmask"))
sras <- read_tsv("data/SraRunTable/SraRunTable.tsv", skip = 1)
sras %>% pull(Organism) %>% unique -> sra.list
genomes.all <- pt %>% filter(Species %in% sra.list) %>% pull(Genome) %>% unique

genomes.all %>% write_lines("data/SRAGenomesList.txt")


genomes.downloaded.fasta <- dir("data/genome", pattern = ".fa$") %>% str_remove(".fa")

genomes.download.2bit <- dir("data/2bit", pattern = ".2bit$") %>% str_remove(".2bit")

genomes.downloaded <- union(genomes.download.2bit, genomes.downloaded.fasta)

genomes.missing <- setdiff(genomes.all, genomes.downloaded)

genomes.missing %>% write_lines("data/missingGenomes.txt")
