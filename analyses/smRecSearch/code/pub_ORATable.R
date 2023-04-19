#!/usr/bin/env Rscript

library(data.table)
library(stringr)
library(magrittr)
source("code/generalFunctions.R")


# ORA.all.f <- "../output/pubFiles/ORA_table_pub.csv"
# ORA.root.dir <- "../output/ORA/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata_RBB_filtered_dyn/pathway_Reactome-top100/"
ORA.all.f <- snakemake@output[[1]]
ORA.root.dir <- snakemake@input[[1]]

dir.create(dirname(ORA.all.f), showWarnings = F, recursive = T)


# ORA.root.dir <- "../output/ORA"
# ORA.root.dir %>% 
#   lineage.pathway.fdr.table %>% 
#   mutate(
#     EnrichRatio = as.numeric(EnrichRatio), FDR=as.numeric(FDR), Size=as.numeric(Size),
#     GeneList = Genes %>% str_replace_all(";", "\n") %>% str_trunc(100, "center", ellipsis = "...\n"),
#     Class = class.pathway(Pathway)
#   ) %>% 
#   fwrite(path=ORA.all.f)


a <- dir(ORA.root.dir, pattern="Project_*", full.names = T) %>% 
  (function(x){x[str_detect(x, "_stable_", negate = T)]}) %>% 
  dir(pattern="enrichment_results", full.names=T) %>% 
  set_names(., str_extract(., "(?<=enrichment_results_(increased|genesDuplicated)_).+(?=.txt)"))
b <- lapply(
  names(a), 
  function(n, df=a){fread(df[[n]])[, branch:=n][]}
) %>% 
  rbindlist
# b %>% split(., .$branch) %>% write.xlsx(ORA.all.f)
b %>% 
  fwrite(ORA.all.f)
