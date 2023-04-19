#!/usr/bin/env Rscript

require(tidyverse)
require(pafr)

get_flip_table <- function(paf){
    rel <- paf %>% 
        group_by(qname, tname) %>% 
        summarize(alen=sum(alen)) %>% 
        group_by(qname) %>% 
        mutate(
               total_alen = sum(alen)
               ) %>% 
        ungroup %>% 
        group_by(qname) %>% 
        mutate(
               rel_alen = alen/total_alen
               ) %>% 
        mutate(
               rel_alen_norm = rel_alen/min(rel_alen)
               ) %>% 
        ungroup() 
    strand <- paf %>%
        group_by(qname, tname) %>% 
        summarize(
                  right=sum(strand=="+"), 
                  wrong = sum(strand=="-")
                  ) 
  rel.top <- rel %>%
      group_by(qname) %>% 
      filter(rel_alen==max(rel_alen)) %>% 
      ungroup 

  full.align <- left_join(rel.top, strand, by=c("qname", "tname")) %>%
     mutate(rel.right = right*rel_alen,
            rel.wrong = rel_alen*wrong, 
            flip = rel.wrong > rel.right) 

  full.align.withnum <- full.align %>%
      mutate(qnum=qname %>% 
                  str_remove("[A-Za-z]+_+") %>% 
                  as.numeric(),
             tnum=tname %>% 
                  str_remove("[A-Za-z]+_+") %>% 
                  as.numeric()
             ) %>% 
      select(qname, qnum, tname, tnum, flip) %>%
      arrange(qnum)
  
  #full.align.sorted <- full.align.withnum %>%
  #    mutate(qname = factor(qname, levels=qnum)) %>% 
  #    arrange(qname) %>% 
  #    select(qname, tname, flip) 
  #return(full.align.sorted)
  return(full.align.withnum)
 }

print(getwd())
print(snakemake@input[[1]])

ali <- read_paf(snakemake@input[[1]])

flip_tab <- ali %>% get_flip_table

flip_tab %>% write_csv(snakemake@output[[1]])

