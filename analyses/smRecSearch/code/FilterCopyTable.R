suppressMessages(suppressWarnings(require(dplyr, quietly = T, warn.conflicts = F)))
suppressMessages(suppressWarnings(require(readr, quietly = T, warn.conflicts = F)))
suppressMessages(suppressWarnings(require(stringr, quietly = T, warn.conflicts = F)))
suppressMessages(suppressWarnings(require(plyranges, quietly = T, warn.conflicts = F)))
suppressMessages(suppressWarnings(source("code/generalFunctions.R")))

##############################################################
genes <- read_lines(snakemake@input[["geneListFiltered"]]) %>%
  #filter_func %>% 
  #pull(Name) %>% 
  unique

print(snakemake@input[["exonTable"]])
read_tsv(snakemake@input[["exonTable"]]) %>% 
  filter(Gene %in% genes) %>% 
  write_tsv(snakemake@output[["exonTable"]])

print(snakemake@input[["copyTable"]])
read_tsv(snakemake@input[["copyTable"]]) %>% 
  filter(Gene %in% genes) %>% 
  write_tsv(snakemake@output[["copyTable"]])

print(snakemake@input[["finalBed"]])
read_tsv(snakemake@input[["finalBed"]]) %>% 
  mutate(Gene=str_remove(name, "_\\d+")) %>% 
  filter(Gene %in% genes) %>% 
  select(-Gene) %>% 
  select(seqnames, start, end, name, score, strand, itemRgb, thickStart, thickEnd, blockCount, blockStart, blockWidth) %>% 
  write_tsv(snakemake@output[["finalBed"]])

print(snakemake@input[["lociTable"]])
read_tsv(snakemake@input[["lociTable"]]) %>% 
  filter(Gene %in% genes) %>% 
  write_tsv(snakemake@output[["lociTable"]])

print(snakemake@input[["duplicateGeneList"]])
read_lines(snakemake@input[["duplicateGeneList"]]) %>% 
  intersect(., genes) %>% 
  write_lines(snakemake@output[["duplicateGeneList"]])
