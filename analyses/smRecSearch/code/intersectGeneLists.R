suppressMessages(suppressWarnings(require(readr)))
suppressMessages(suppressWarnings(require(dplyr)))

a <- read_lines(snakemake@input[["geneListA"]])
b <- read_lines(snakemake@input[["geneListB"]])

setdiff(a, b) %>% 
  write_lines(snakemake@output[["A_not_B"]])
setdiff(b, a) %>% 
  write_lines(snakemake@output[["B_not_A"]])
intersect(a, b) %>% 
  write_lines(snakemake@output[["A_and_B"]])
