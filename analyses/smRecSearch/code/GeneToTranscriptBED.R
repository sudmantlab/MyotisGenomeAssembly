require(dplyr, quietly = T)
require(tidyr, quietly = T)
require(magrittr, quietly = T)
require(readr, quietly = T)
require(plyranges, quietly = T)
source("../code/generalFunctions.R")

##############################################################

print("Reading Translation Table...")
translation.table <- read_tsv(snakemake@input[["TranslationTable"]])
print("Done.")

print("Reading RBB file...")
rbb <- get.rbb(snakemake@input[["RBB_file"]])%>%
  mutate(Gene = sapply(QName, function(id, tt = translation.table){tt$Name[which(tt$ID == id)][1]}))
print("Done.")

print("Getting BED entries for gene of interest...")
rbb.subset <- rbb %>% 
  filter(Gene == snakemake@wildcards[["gene"]]) %>% 
  mutate(name = paste(Gene, name, sep="."))
print("Done.")

print("Writing BED file...")
rbb.subset %>% 
  as_tibble %>% 
  write_bed(snakemake@output[1])