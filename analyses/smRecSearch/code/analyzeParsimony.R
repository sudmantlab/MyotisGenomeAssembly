#suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(ape)))
suppressMessages(suppressWarnings(require(treeio)))
suppressMessages(suppressWarnings(require(tidytree)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(readr)))
suppressMessages(suppressWarnings(require(magrittr)))
suppressMessages(suppressWarnings(require(stringr)))
source("code/generalFunctions.R")

dir.create("output/geneLists/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/parsimony-strict", showWarnings = F)


# Note: these are for strict parsimony

has.increased <- function(gene, a, anc) {
  copy.anc.high <- anc[[paste0(gene, ".max")]] 
  copy.anc.low <- anc[[paste0(gene, ".max")]] 
  copy.a <- a[[paste0(gene, ".min")]] 
  if (is.na(copy.a) | is.na(copy.anc.high)| is.na(copy.anc.low)){
    return(FALSE)
  } else{
    return(copy.a > copy.anc.high)
  }
}

has.stable <- function(gene, a, anc) {
  copy.anc.high <- anc[[paste0(gene, ".max")]] 
  copy.anc.low <- anc[[paste0(gene, ".max")]] 
  copy.a <- a[[paste0(gene, ".min")]] 
  if (is.na(copy.a) | is.na(copy.anc.high)| is.na(copy.anc.low)){
    return(FALSE)
  } else{
    return(copy.anc.low <= copy.a & copy.a <= copy.anc.high)
  }
}

has.decreased <- function(gene, a, anc) {
  copy.anc.high <- anc[[paste0(gene, ".max")]] 
  copy.anc.low <- anc[[paste0(gene, ".max")]] 
  copy.a <- a[[paste0(gene, ".min")]] 
  if (is.na(copy.a) | is.na(copy.anc.high)| is.na(copy.anc.low)){
    return(FALSE)
  } else{
    return(copy.a < copy.anc.low)
  }
}


run.comparison <- function(comparison, type, rootdir, l=afrotheria.list, gl=genelist){
  ancest <- l[[comparison[1]]]
  terminal <- l[[comparison[2]]]
  if (type=="increased"){
    v <- gl %>% 
      lapply(
        ., 
        . %>% has.increased(., a=terminal, anc=ancest)
      ) %>% 
      .[which(.==T)] %>% 
      names
  } else if (type == "decreased"){
    v <- gl %>% 
      lapply(
        ., 
        . %>% has.decreased(., a=terminal, anc=ancest)
      ) %>% 
      .[which(.==T)] %>% 
      names
  } else if (type == "stable"){
    v <- gl %>% 
      lapply(
        ., 
        . %>% has.stable(., a=terminal, anc=ancest)
      ) %>% 
      .[which(.==T)] %>% 
      names
  } else{
    return(NA)
  }
  f <- paste0(rootdir, type, "-", comparison[1], "-to-", comparison[2], ".txt")
  write_lines(v, f)
}

coalesce_join <- function(x, y, 
                          by = "node", suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  # Source: https://alistaire.rbind.io/blog/coalescing-joins/
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  #printCoalesce <- function(x, s=suffix){
  #  print(x)
  #  dplyr::coalesce(
  #    joined[[paste0(x, suffix[1])]],
  #    joined[[paste0(x, suffix[2])]]
  #  )}
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    as.double(joined[[paste0(.x, suffix[1])]]),
    as.double(joined[[paste0(.x, suffix[2])]])
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}


parsimony.treedata <- read.beast("output/parsimony/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata_RBB_filtered-parsimony.nexus")
# parsimony.treedata <- read.beast(snakemake@input[[1]])
tree <- parsimony.treedata@phylo

genelist <- names(parsimony.treedata@data) %>% str_remove_all(".min|.max|node") %>% unique %>% stringi::stri_remove_empty() %>% set_names(.,.)

# test.gene <- genelist[1]

parsimony.treedata <- parsimony.treedata %>% as_tibble
#parsimony.treedata <- parsimony.treedata %>% as.treedata()

parsimony.treedata[which(parsimony.treedata$label==""),] <- "Root"

afrotheria.list <- list(
  "Loxodonta.africana"= parsimony.treedata %>% filter(label=="Loxodonta.africana"),
  "Trichechus.manatus"=parsimony.treedata %>% filter(label=="Trichechus.manatus"),
  "Procavia.capensis"= parsimony.treedata %>% filter(label=="Procavia.capensis"),
  "Orycteropus.afer"= parsimony.treedata %>% filter(label=="Orycteropus.afer"),
  "Dasypus.novemcinctus"= parsimony.treedata %>% filter(label=="Choloepus.hoffmanni"),
  "Choloepus.hoffmanni" = parsimony.treedata %>% filter(label=="Dasypus.novemcinctus"),
  "Chrysochloris.asiatica" = parsimony.treedata %>% filter(label=="Chrysochloris.asiatica"),
  "Echinops.telfairi" = parsimony.treedata %>% filter(label=="Echinops.telfairi"),
  "Elephantulus.edwardii" = parsimony.treedata %>% filter(label=="Elephantulus.edwardii"),
  "Loxodonta.cyclotis" = parsimony.treedata %>% filter(label=="Loxodonta.cyclotis"),
  "Palaeoloxodon.antiquus" = parsimony.treedata %>% filter(label=="Palaeoloxodon.antiquus"),
  "Elephas.maximus" = parsimony.treedata %>% filter(label=="Elephas.maximus"),
  "Mammuthus.columbi" = parsimony.treedata %>% filter(label=="Mammuthus.columbi"),
  "Mammuthus.primigenius" = parsimony.treedata %>% filter(label=="Mammuthus.primigenius"),
  "Mammut.americanum" = parsimony.treedata %>% filter(label=="Mammut.americanum"),
  "Xenarthra" = parsimony.treedata %>% filter(label=="Xenarthra"),
  "Afrotheria" =  parsimony.treedata %>% filter(label=="Afrotheria"),
  "Afroinsectivora" =  parsimony.treedata %>% filter(label=="Afroinsectivora"),
  "Afrosoricida" =  parsimony.treedata %>% filter(label=="Afrosoricida"),
  "Pseudoungulata" =  parsimony.treedata %>% filter(label=="Pseudoungulata"),
  "Paenungulata" =  parsimony.treedata %>% filter(label=="Paenungulata"),
  "Tehytheria" =  parsimony.treedata %>% filter(label=="Tehytheria"),
  "Proboscidea" =  parsimony.treedata %>% filter(label=="Proboscidea"),
  "Elephantidae" =  parsimony.treedata %>% filter(label=="Elephantidae"),
  "Loxodontini" =  parsimony.treedata %>% filter(label=="Loxodontini"),
  "Loxodona" =  parsimony.treedata %>% filter(label=="Loxodona"),
  "Elephantina" =  parsimony.treedata %>% filter(label=="Elephantina"),
  "Mammuthus" =  parsimony.treedata %>% filter(label=="Mammuthus"),
  "Root" = parsimony.treedata %>% filter(label=="Root")
)

comparisons <- list(
  c("Proboscidea", "Elephantidae"),
  c("Elephantidae", "Elephantina"),
  c("Elephantina", "Elephas.maximus"),
  c("Elephantina", "Mammuthus"),
  c("Mammuthus", "Mammuthus.columbi"),
  c("Mammuthus", "Mammuthus.primigenius"),
  c("Elephantidae", "Loxodontini"),
  c("Loxodontini", "Loxodonta.africana"),
  c("Loxodontini", "Loxodona"),
  c("Loxodona", "Loxodonta.cyclotis"),
  c("Loxodona", "Palaeoloxodon.antiquus"),
  c("Proboscidea", "Mammut.americanum"),
  c("Tehytheria", "Proboscidea"),
  c("Tehytheria", "Trichechus.manatus"),
  c("Paenungulata", "Tehytheria"),
  c("Paenungulata", "Procavia.capensis"),
  c("Pseudoungulata", "Paenungulata"),
  c("Pseudoungulata", "Orycteropus.afer"),
  c("Afrotheria", "Pseudoungulata"),
  c("Afrotheria", "Afroinsectivora"),
  c("Afroinsectivora", "Elephantulus.edwardii"),
  c("Afroinsectivora", "Afrosoricida"),
  c("Afrosoricida", "Echinops.telfairi"),
  c("Afrosoricida", "Chrysochloris.asiatica"),
  c("Root", "Afrotheria"),
  c("Root", "Choloepus.hoffmanni"),
  c("Root", "Dasypus.novemcinctus")
) %>% set_names(., sapply(.,. %>% paste(., collapse="-to-")))
root.dir <- "output/geneLists/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/parsimony-strict/"
for (compare in comparisons){
  run.comparison(comparison = compare, type = "increased", l=afrotheria.list, gl=genelist, rootdir = root.dir)
  # run.comparison(comparison = compare, type = "decreased", l=afrotheria.list, gl=genelist, rootdir = root.dir)
  run.comparison(comparison = compare, type = "stable", l=afrotheria.list, gl=genelist, rootdir = root.dir)
}

