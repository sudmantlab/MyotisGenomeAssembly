suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(magrittr)))
suppressMessages(suppressWarnings(require(stringr)))
suppressMessages(suppressWarnings(require(readr)))
suppressMessages(suppressWarnings(require(UpSetR)))
source("code/generalFunctions.R")

is.cancer <- function(row) {
  str_detect(row[["Gene"]], "APC|[Cc]yclin|[Cc]ell[Cc]ycle|Ubiq|NER|ROS|DNA|G1|[Mm]eiotic|[Mm]itotic|[Cc]heckpoint|[Cc]ancer")
}

read_tsv_if <- function(x){
  if(length(x) == 0){
    tibble(geneSet=character(0), description=character(0), link=character(0), size=character(0), overlap=character(0), 
           expect=character(0), enrichmentRatio=character(0), pValue=character(0), FDR=character(0), overlapId=character(0), 
           userId=character(0))
  } else{
    read_tsv(x) 
  }
}

comparisons.ml <- dir("output/ORA/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/") %>% set_names(.,.)
tables.ml <- comparisons %>% 
  lapply(
  ., 
  . %>% 
    paste0("output/ORA/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/", .) %>% 
    dir(., full.names = T) %>% 
    dir(., full.names = T, pattern = "enrichment_results.*") %>% 
    read_tsv_if
)
comparisons.pr <- dir("output/ORA/parsimony-strict/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/") %>% set_names(.,.)
tables.pr <- comparisons %>% 
  lapply(
    ., 
    . %>% 
      paste0("output/ORA/parsimony-strict/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/", .) %>% 
      dir(., full.names = T) %>% 
      dir(., full.names = T, pattern = "enrichment_results.*") %>% 
      read_tsv_if
  )

bin.ml <- tables.ml %>% 
  lapply(
    .,
    . %>% pull(description)
  ) %>% 
  binTableFromList() %>% 
  mutate(
    cancer=str_detect(Gene,"APC|[Cc]yclin|[Cc]ell[Cc]ycle|Ubiq|NER|ROS|DNA|G1|[Mm]eiotic|[Mm]itotic|[Cc]heckpoint|[Cc]ancer|E[123]")
    )

bin.pr <- tables.pr %>% 
  lapply(
    .,
    . %>% pull(description)
  ) %>% 
  binTableFromList()

names(bin.ml %>% select(-Gene)) %>% sort

decreased <- names(bin.ml) %>% str_subset(pattern = "decreased")
# decreased %>% sapply(., . %>% paste0("list(\nquery=intersects,\nparams=c('", ., "'),\ncolor='red',\nactive=T),\n")) %>% cat()
increased <- names(bin.ml) %>% str_subset(pattern = "increased")
stable <- names(bin.ml) %>% str_subset(pattern = "stable")

upset(bin.ml, keep.order = T, 
      queries = list(
        list(
          query=intersects,
          params=c('decreased-Paenungulata-to-Procavia.capensis'),
          color='red',
          active=T),
        list(
          query=intersects,
          params=c('decreased-Pseudoungulata-to-Loxodonta'),
          color='red',
          active=T),
        list(
          query=intersects,
          params=c('decreased-Pseudoungulata-to-Loxodonta.africana'),
          color='red',
          active=T)
      )
      )

      