library(tidyverse)

f.diamond <- snakemake@input[[1]]
# f.diamond <- 'output/EVM_manual/mMyoAui1_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/mMyoAui1_TOGARegOnly_Augustus_MiniprotMerged_GenBankLiftoff.complementedTOGA.diamond.tsv'

o.diamond.crit <- snakemake@output[['tbl']]
# o.diamond.crit <- 'output/EVM_manual/mMyoAui1_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/mMyoAui1_TOGARegOnly_Augustus_MiniprotMerged_GenBankLiftoff.complementedTOGA.diamond.criteria.tsv'

o.diamond.keepers <- snakemake@output[['keep']]
# o.diamond.keepers <- 'output/EVM_manual/mMyoAui1_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/mMyoAui1_TOGARegOnly_Augustus_MiniprotMerged_GenBankLiftoff.complementedTOGA.diamond.keep.ids'


diamond <- read_tsv(
  f.diamond, 
  col_names = c('qseqid', 'staxids', 'bitscore', 'sseqid', 
                'pident', 'length', 'mismatch', 'gapopen', 
                'qlen', 'qstart', 'qend', 'slen', 'sstart', 
                'send', 'ppos', 'evalue')
  ) %>% 
  mutate(qcov = (qend-qstart)/qlen, scov = (send-sstart)/slen)

diamond.crit <- diamond %>% 
  transmute(
    qseqid = qseqid, 
    sseqid = sseqid, 
    'qcov>=75' = (qcov >= 0.75), 
    'scov>=50' = (scov >= 0.5), 
    'ppos>=50' = (ppos>=50.0)
  )

diamond.crit %>% write_tsv(o.diamond.crit)

diamond.keepers <- diamond.crit %>% 
  filter(
    `qcov>=75`, 
    `scov>=50`, 
    `ppos>=50`
  ) %>% 
  pull(qseqid)

diamond.keepers %>% write_lines(o.diamond.keepers)
