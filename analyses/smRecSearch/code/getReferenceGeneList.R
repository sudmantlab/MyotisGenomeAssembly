suppressMessages(suppressWarnings(require(tidyverse)))
source("code/generalFunctions.R")

translation.table <- read_tsv(snakemake@input[["TranslationTable"]])

translation.table %>%
  pull(Name) %>%
  unique %>%
  write_lines(snakemake@output[["geneList"]])

translation.table %>%
  filter_func %>%
  pull(Name) %>%
  unique %>%
  write_lines(snakemake@output[["filteredGeneList"]])
