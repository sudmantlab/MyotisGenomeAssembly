suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(ape)))
suppressMessages(suppressWarnings(library(treeio)))
suppressMessages(suppressWarnings(library(tidytree)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(tidyr)))
#source("code/generalFunctions.R")
# library(ggtree)

dynToCN <- function(val){
  val <- as.character(val)
  scale = c(NA, 1:36) %>% 
    set_names(c("-", as.character(0:9), LETTERS))
  if ( ! val %in% names(scale)) {stop(sprintf("%s not in scale!", val))}
  newval <- scale[[val]]
  #print(val)
  return(newval)
}

dynToCN2 <- function(val){
  val <- as.character(val)
  scale = c(NA, 1:36) %>% 
    set_names(c("?", as.character(0:9), LETTERS))
  if ( ! val %in% names(scale)) {stop(sprintf("%s not in scale!", val))}
  newval <- scale[[val]]
  # print(length(val))
  # print(class(val))
  # print(val)
  return(as.numeric(scale[[val]]))
}


# Note: these are for ML

has.increased <- function(gene, a, anc) {
  copy.anc <- anc[[gene]] 
  copy.a <- a[[gene]] 
  if (is.na(copy.a) | is.na(copy.anc)){
    return(FALSE)
  } else{
    return(copy.a > copy.anc)
  }
}

has.stable <- function(gene, a, anc) {
  copy.anc <- anc[[gene]] 
  copy.a <- a[[gene]] 
  if (is.na(copy.a) | is.na(copy.anc)){
    return(FALSE)
  } else{
    return(copy.a == copy.anc)
  }
}

has.decreased <- function(gene, a, anc) {
  copy.anc <- anc[[gene]] 
  copy.a <- a[[gene]] 
  if (is.na(copy.a) | is.na(copy.anc)){
    return(FALSE)
  } else{
    return(copy.a < copy.anc)
  }
}


run.comparison <- function(comparison, type, rootdir, l, gl){
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


dir.create(snakemake@params[["dir"]], showWarnings = F, recursive = T)

tree.genechanges <- read.tree(snakemake@input[["iqtree"]])%>% 
  root("Dasypus.novemcinctus")
# tree.genechanges <- read.tree("output/iqtree/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/atlantogenata_RBB_filtered_dyn.treefile") %>%
#   root("Dasypus.novemcinctus")

# tree.genechanges %>% ggtree() + geom_tiplab() + theme_tree2() + xlim(0,.6)

tree.states <- read_tsv(snakemake@input[["states"]], comment = "#") %>%
  mutate(State = State %>% sapply(dynToCN))
# tree.states <- read_tsv("output/iqtree/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/atlantogenata_RBB_filtered_dyn.state", comment = "#") %>%
#   mutate(State = State %>% sapply(dynToCN))
# tree.states

geneCopyTable.wide <- read_csv(snakemake@input[["wideTable_dyn"]]) %>%
  group_by(label) %>% 
  mutate_at(vars(-group_cols()), . %>% sapply(dynToCN2)) %>% 
  ungroup
# geneCopyTable.wide <- read_csv("output/geneCopyTable/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata-GeneCopyTable_RBB_filtered_wide_dyn.csv") %>%
#   group_by(label) %>%
#   mutate_at(vars(-group_cols()), . %>% sapply(dynToCN2)) %>%
#   ungroup

sitenames <- read_lines(snakemake@input[["sitenames"]])%>% 
  tibble(Gene = ., Site = seq_along(.))
# sitenames <- read_lines("output/phylip/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata_hq-GeneCopyTable_RBB_filtered_coded_list.txt") %>%
#   tibble(Gene = ., Site = seq_along(.))

# length(sitenames$Site) == length(unique(tree.states$Site))

tree.states %<>% full_join(., sitenames, by="Site") %>% select(-Site) %>% rename(label=Node)
# tree.states %>% pull(node) %>% unique
# tree.genechanges$node.label %>% unique

tree.states.wide <- tree.states %>% 
  select(label, Gene, State) %>% 
  pivot_wider(names_from = Gene, values_from = State)

tree.states.wide.full <- rbind(geneCopyTable.wide, tree.states.wide) %>% 
  #mutate_if(is.character, list(~replace_na(., "-")))
  mutate_if(is.character, list(~na_if(., "-")))

treedata.genechanges <- full_join(tree.genechanges, tree.states.wide.full, by="label")


# 
# treedata.genechanges %>% 
#   ggtree() + 
#   geom_tiplab() + 
#   geom_nodelab()

write.beast(treedata.genechanges, snakemake@output[["nexus"]])
# write.beast(treedata.genechanges, "output/geneLists/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata_RBB_filtered_dyn/atlantogenata_RBB_filtered_dyn.nexus")


tree.ml <- treedata.genechanges@phylo
maxLike.treedata <- treedata.genechanges %>% as_tibble
#maxLike.treedata <- maxLike.treedata %>% as.treedata()

afrotheria.list.ml <- list(
  "Loxodonta.africana"= maxLike.treedata %>% filter(label=="Loxodonta.africana"),
  "Trichechus.manatus"=maxLike.treedata %>% filter(label=="Trichechus.manatus"),
  "Procavia.capensis"= maxLike.treedata %>% filter(label=="Procavia.capensis"),
  "Orycteropus.afer"= maxLike.treedata %>% filter(label=="Orycteropus.afer"),
  "Dasypus.novemcinctus"= maxLike.treedata %>% filter(label=="Choloepus.hoffmanni"),
  "Choloepus.hoffmanni" = maxLike.treedata %>% filter(label=="Dasypus.novemcinctus"),
  "Chrysochloris.asiatica" = maxLike.treedata %>% filter(label=="Chrysochloris.asiatica"),
  "Echinops.telfairi" = maxLike.treedata %>% filter(label=="Echinops.telfairi"),
  "Elephantulus.edwardii" = maxLike.treedata %>% filter(label=="Elephantulus.edwardii"),
  "Loxodonta.cyclotis" = maxLike.treedata %>% filter(label=="Loxodonta.cyclotis"),
  "Palaeoloxodon.antiquus" = maxLike.treedata %>% filter(label=="Palaeoloxodon.antiquus"),
  "Elephas.maximus" = maxLike.treedata %>% filter(label=="Elephas.maximus"),
  "Mammuthus.columbi" = maxLike.treedata %>% filter(label=="Mammuthus.columbi"),
  "Mammuthus.primigenius" = maxLike.treedata %>% filter(label=="Mammuthus.primigenius"),
  "Mammut.americanum" = maxLike.treedata %>% filter(label=="Mammut.americanum"),
  "Xenarthra" = maxLike.treedata %>% filter(label=="Xenarthra"),
  "Afrotheria" =  maxLike.treedata %>% filter(label=="Afrotheria"),
  "Afroinsectivora" =  maxLike.treedata %>% filter(label=="Afroinsectivora"),
  "Afrosoricida" =  maxLike.treedata %>% filter(label=="Afrosoricida"),
  "Pseudoungulata" =  maxLike.treedata %>% filter(label=="Pseudoungulata"),
  "Paenungulata" =  maxLike.treedata %>% filter(label=="Paenungulata"),
  "Tethytheria" =  maxLike.treedata %>% filter(label=="Tethytheria"),
  "Proboscidea" =  maxLike.treedata %>% filter(label=="Proboscidea"),
  "Elephantidae" =  maxLike.treedata %>% filter(label=="Elephantidae"),
  "Loxodontini" =  maxLike.treedata %>% filter(label=="Loxodontini"),
  "Loxodona" =  maxLike.treedata %>% filter(label=="Loxodona"),
  "Elephantina" =  maxLike.treedata %>% filter(label=="Elephantina"),
  "Mammuthus" =  maxLike.treedata %>% filter(label=="Mammuthus"),
  "Root" = maxLike.treedata %>% filter(label=="Root")
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
  c("Tethytheria", "Proboscidea"),
  c("Tethytheria", "Trichechus.manatus"),
  c("Paenungulata", "Tethytheria"),
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

genelist <- maxLike.treedata %>% select(-label, -node,-branch.length,-parent) %>% names %>% set_names(.,.)

for (compare in comparisons){
  print(compare)
  run.comparison(comparison = compare, type = "increased", l=afrotheria.list.ml, gl=genelist, rootdir = snakemake@params[["dir"]])
  # run.comparison(comparison = compare, type = "decreased", l=afrotheria.list, gl=genelist, rootdir = root.dir)
  run.comparison(comparison = compare, type = "stable", l=afrotheria.list.ml, gl=genelist, rootdir = snakemake@params[["dir"]])
}
