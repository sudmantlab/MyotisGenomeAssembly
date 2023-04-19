suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(ape)))
suppressMessages(suppressWarnings(require(treeio)))
suppressMessages(suppressWarnings(require(tidytree)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(readr)))
suppressMessages(suppressWarnings(require(magrittr)))
suppressMessages(suppressWarnings(require(stringr)))
suppressMessages(suppressWarnings(require(tidyr)))
source("code/generalFunctions.R")
# library(ggtree)


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


dir.create("output/geneLists/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/", showWarnings = F)

# tree.genechanges <- read.tree(snakemake@input[["rateTree"]])
tree.genechanges <- read.tree("output/iqtree/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/atlantogenata_evidenced_filtered.treefile") %>% 
  root("Dasypus.novemcinctus")

# tree.genechanges %>% ggtree() + geom_tiplab() + theme_tree2() + xlim(0,.6)

# tree.states <- read_table(snakemake@input[["state"]])
tree.states <- read_tsv("output/iqtree/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/atlantogenata_evidenced_filtered.state", comment = "#")
# tree.states

# geneCopyTable.wide <- fread(snakemake@input[["wideTable"]])
geneCopyTable.wide <- read_csv("output/geneCopyTable/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata-GeneCopyTable_evidenced_filtered_wide.csv")

# sitenames <- read_lines(snakemake@input[["siteNames"]])
sitenames <- read_lines("output/phylip/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata-GeneCopyTable_evidenced_filtered_coded_list.txt") %>% 
  tibble(Gene = ., Site = seq_along(.))

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

write.beast(treedata.genechanges, "output/geneLists/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/atlantogenata_evidenced_filtered.nexus")

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
  "Tehytheria" =  maxLike.treedata %>% filter(label=="Tehytheria"),
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

genelist <- maxLike.treedata %>% select(-label, -node,-branch.length,-parent) %>% names %>% set_names(.,.)

root.dir.ml <- "output/geneLists/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/"
for (compare in comparisons[19:length(comparisons)]){
  print(compare)
  run.comparison(comparison = compare, type = "increased", l=afrotheria.list.ml, gl=genelist, rootdir = root.dir.ml)
  # run.comparison(comparison = compare, type = "decreased", l=afrotheria.list, gl=genelist, rootdir = root.dir)
  run.comparison(comparison = compare, type = "stable", l=afrotheria.list.ml, gl=genelist, rootdir = root.dir.ml)
}
