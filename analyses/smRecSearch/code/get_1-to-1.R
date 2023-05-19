library(tidyverse)

df <- read_tsv(snakemake@input[[1]],
               col_names = c("chr", "s", "e", "id_full", "score", "strand", "tstart", "tend", "color", "nblock", "block_len", "block_start", 
                             "rhit", "rstart", "rend", "qcov"),
               col_types = c("chr"="c", "s"="n", "e"="n", "id_full"="c", "score"='n', "strand"="c", "tstart"="c", "tend"="c", "color"="c", 
                             "nblock"="c", "block_len"="c", "block_start"="c", "rhit"="c", "rstart"="c", "rend"="c", "qcov"="c")) %>%
  select(-rhit, -rstart, -rend, -qcov) %>% 
  arrange_all() %>% 
  distinct 

main_hits <- df %>% 
  select(id_full, score) %>% 
  separate(id_full, c("id", "copy"), sep="_") %>% 
  group_by(id) %>%
  filter(score==max(score)) %>%
  filter(copy==min(copy)) %>%
  ungroup %>% 
  unite("id_full", id, copy, sep="_") %>% 
  group_by(id_full) %>%
  summarize(id_full = first(id_full)) %>%
  ungroup %>%
  pull(id_full) %>%
  unique

df.final <- df %>% 
  filter(id_full %in% main_hits) %>% 
  group_by(id_full) %>%
  filter(row_number()==1)

df.final %>% write_tsv(snakemake@output[[1]], col_names=F)
