library(tidyverse)

f.diamond <- snakemake@input[['diamond']]

f.gff <- snakemake@input[['gff']]

f.out.tsv <- snakemake@output[['tsv']]

f.out.txt <- snakemake@output[['txt']]

diamond <- read_tsv(
  f.diamond,
  col_names = c('qseqid', 'bitscore', 'sseqid',
                'pident', 'length', 'mismatch', 'gapopen',
                'qlen', 'qstart', 'qend', 'slen', 
                'sstart', 'send', 'ppos', 'evalue')
  )

genes <- read_tsv(f.gff, col_names=F, comment="#") %>% 
  filter(X3 == "mRNA") %>% 
  transmute(
    qseqid = X9 %>% 
               str_extract("ID=[A-Za-z0-9]+;") %>% 
               str_remove("ID=") %>% 
               str_remove(";"), 
    Target = X9 %>% 
               str_extract("Target=[-_A-Za-z0-9.]+ ") %>% 
               str_remove("Target=") %>% 
               str_remove(" ")
    )

diamond.full <- diamond %>%
  left_join(genes, by='qseqid')

diamond.full %>% select(qseqid, Target, sseqid) %>% head %>% print

diamond.rbh <- diamond.full %>%
  filter(sseqid == Target)

diamond.rbh %>% nrow %>% print

diamond.list <- diamond.rbh %>%
  pull(qseqid) %>%
  unique()

#nrow(diamond.rbh) %>% print
diamond.rbh %>% write_tsv(f.out.tsv)
diamond.list %>% write_lines(f.out.txt)
