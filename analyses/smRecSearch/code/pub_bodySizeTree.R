#!/usr/bin/env Rscript

library(tidyverse)
library(treeio)
library(tidytree)
library(ggtree)

# bsize.tree.f <- "../output/AncSizeRec/AnnotatedTreeForFigure.tree"
# nexus.f <- "../data/stableTraits/eutheria.nexus"
# ancstates.f <- "../data/stableTraits/eutheria.ancstates"
# brlens.f <- "../data/stableTraits/eutheria.brlens"
# genomeTable.f <- "../output/other/genomeTable.csv"

bsize.tree.f <- snakemake@input[["bsize"]]
nexus.f <- snakemake@input[["nexus"]]
ancstates.f <- snakemake@input[["ancstates"]]
brlens.f <- snakemake@input[["brlens"]]
genomeTable.f <- snakemake@input[["genomeTable"]]
  
# outTree.f <- "../output/pubFiles/bodysize-Eutheria.csv"
# outTree.Atlantogenata.f <- "../output/pubFiles/bodysize-Atlantogenata.csv"

outTree.f <- snakemake@output[["outTree_full"]]
outTree.Atlantogenata.f <- snakemake@output[["outTree_Atlantogenata"]]

dir.create(dirname(outTree.f), showWarnings = F, recursive = T)

genomeTable <- read_csv(genomeTable.f)

bsize.tree <- read.beast(bsize.tree.f) %>% 
  as_tibble %>% 
  mutate(sqrt.rate = branch.length) %>% 
  as.treedata
# st.sqrt <- read.tree("../data/stableTraits/eutheria.sqchange.tree") 
# st.sqrt$node.label <- st.sqrt$node.label %>% str_remove("node")
# st.sqrt %<>% 
#   as_tibble %>% 
#   mutate(sqrt.rate = branch.length) %>% 
#   as.treedata
st.nexus <- read.nexus(nexus.f) %>%
  as_tibble() %>%
  rename(real_label=label) %>%
  mutate(
    name = str_extract(real_label, "[A-Za-z][A-Za-z0-9*.]+'?$"),
    label = str_extract(real_label, "^(\\d)+(?=\\.[:alpha:])|^'(\\d)+(?=\\.[:alpha:])|^(\\d)+$") %>% str_remove("'")
  )
st.nexus[is.na(st.nexus$label),]$label <- st.nexus[is.na(st.nexus$label),]$name
classes.nexus <- st.nexus %>% class
st.nexus %<>% as.treedata()
st.other <- full_join(
  read_tsv(ancstates.f) %>% mutate(Parameter = str_remove(Parameter, "^n(?=\\d)")) %>% rename(label= Parameter), 
  read_tsv(brlens.f) %>% mutate(label=Branch %>% str_split("->") %>% sapply("[",1) %>% str_remove("^n(?=\\d)")), 
  by="label"
) %>% 
  mutate(node = Branch %>% str_split("->") %>% sapply("[",1) %>% str_remove("^n(?=\\d)"), parent = Branch %>% str_split("->") %>% sapply("[",2) %>% str_remove("^n(?=\\d)"))
st.other[st.other$label=="0",]$node <- "0"
#class(st.other) <- classes.nexus#c(class(st.other), "tbl_tree")
st.other %<>% select(parent, node, label, everything()) %>% filter(!is.na(parent)) %>% arrange(parent, node) %>% as.treedata()
st.part <- merge_tree(st.nexus, st.other)
st.full <- merge_tree(st.part, bsize.tree)

# Tree
st.innerNodes <- tribble(
  ~"Node",                                                     ~"label",
  MRCA(st.full, "Loxodonta.africana", "Loxodonta.cyclotis"), "Loxodontini",
  MRCA(st.full, "Palaeoloxodon.antiquus", "Loxodonta.cyclotis"), "Loxodona",
  MRCA(st.full, "Loxodonta.africana", "Elephas.maximus"), "Elephantidae",
  MRCA(st.full, "Elephas.maximus", "Mammuthus.primigenius"), "Elephantina",
  MRCA(st.full, "Mammuthus.columbi", "Mammuthus.primigenius"), "Mammuthus",
  MRCA(st.full, "Loxodonta.africana", "Mammut.americanum"), "Elephantimorpha",
  MRCA(st.full, "Loxodonta.africana", "Numidotherium.koholense"), "Proboscidea",
  MRCA(st.full, "Loxodonta.africana", "Trichechus.manatus"), "Tethytheria",
  MRCA(st.full, "Loxodonta.africana", "Procavia.capensis"), "Paenungulata",
  MRCA(st.full, "Loxodonta.africana", "Orycteropus.afer"), "Pseudoungulata",
  MRCA(st.full, "Loxodonta.africana", "Echinops.telfairi"), "Afrotheria",
  MRCA(st.full, "Chrysochloris.asiatica", "Echinops.telfairi"), "Afrosoricida",
  MRCA(st.full, "Elephantulus.edwardii", "Chambius.kasserinensis"), "Macroscelidae",
  MRCA(st.full, "Dasypus.novemcinctus", "Choloepus.hoffmanni"), "Xenarthra",
  MRCA(st.full, "Dasypus.novemcinctus", "Loxodonta.africana"), "Atlantogenata"
)
st.full %<>% 
  as_tibble() %>% 
  select(-label.y, label=label.x)

st.full[which(st.full$node %in% st.innerNodes$Node),]$label <- sapply(st.full[which(st.full$node %in% st.innerNodes$Node),]$node, function(n, tb=st.innerNodes){tb %>% filter(Node==n) %>% pull(label) %>% unique})
# A bit of manual labor
st.full[which(st.full$label %in% as.character(1859:1876)),]$label <- NA

st.ancAtlantogenata <- st.innerNodes %>% filter(label=="Atlantogenata") %>% pull(Node)

st.tree <- full_join(st.full %>% as.treedata(), genomeTable, by = "label") %>% as_tibble() %>% filter(!is.na(node)) %>% as.treedata()

st.Atlantogenata <- st.tree %>% 
  # as.treedata() %>% 
  tree_subset(., node=st.ancAtlantogenata, levels_back=0) %>%
  as_tibble() %>% 
  select(parent, node, branch.length, label, Rate, Orig_len, Brownian, Median, `95%CI_Low`, `95%CI_High`, sqrt.rate, lnSize, Genome, Species.strict, `Common Name`, Species, Spc.short, BestGenome, phyloPicUID) %>%
  mutate(branch.length=Orig_len) %>% 
  mutate_at(vars(-label, -Genome, -Species.strict, -`Common Name`, -Species, -Spc.short, -BestGenome, -phyloPicUID), . %>% as.numeric %>% round(2)) %>% 
  mutate(label=label %>% str_replace_all("\\.", " ")) %>% 
  as.treedata()

# st.tree %>% as_tibble() %>% filter(node==ancNode.Atlantogenata)

tree.Atlantogenata <- st.Atlantogenata

tree.Atlantogenata.table <- tree.Atlantogenata %>% as_tibble

tree.Atlantogenata.nodes <- list(
  "Loxodonta.africana"= tree.Atlantogenata.table %>% filter(label=="Loxodonta africana") %>% pull(node) %>% unique,
  "Trichechus.manatus"= tree.Atlantogenata.table %>% filter(label=="Trichechus manatus") %>% pull(node) %>% unique,
  "Procavia.capensis"= tree.Atlantogenata.table %>% filter(label=="Procavia capensis") %>% pull(node) %>% unique,
  "Orycteropus.afer"= tree.Atlantogenata.table %>% filter(label=="Orycteropus afer") %>% pull(node) %>% unique,
  "Dasypus.novemcinctus"= tree.Atlantogenata.table %>% filter(label=="Choloepus hoffmanni") %>% pull(node) %>% unique,
  "Choloepus.hoffmanni" = tree.Atlantogenata.table %>% filter(label=="Dasypus novemcinctus") %>% pull(node) %>% unique,
  "Chrysochloris.asiatica" = tree.Atlantogenata.table %>% filter(label=="Chrysochloris asiatica") %>% pull(node) %>% unique,
  "Echinops.telfairi" = tree.Atlantogenata.table %>% filter(label=="Echinops telfairi") %>% pull(node) %>% unique,
  "Elephantulus.edwardii" = tree.Atlantogenata.table %>% filter(label=="Elephantulus edwardii") %>% pull(node) %>% unique,
  "Loxodonta.cyclotis" = tree.Atlantogenata.table %>% filter(label=="Loxodonta cyclotis") %>% pull(node) %>% unique,
  "Palaeoloxodon.antiquus" = tree.Atlantogenata.table %>% filter(label=="Palaeoloxodon antiquus") %>% pull(node) %>% unique,
  "Elephas.maximus" = tree.Atlantogenata.table %>% filter(label=="Elephas maximus") %>% pull(node) %>% unique,
  "Mammuthus.columbi" = tree.Atlantogenata.table %>% filter(label=="Mammuthus columbi") %>% pull(node) %>% unique,
  "Mammuthus.primigenius" = tree.Atlantogenata.table %>% filter(label=="Mammuthus primigenius") %>% pull(node) %>% unique,
  "Mammut.americanum" = tree.Atlantogenata.table %>% filter(label=="Mammut americanum") %>% pull(node) %>% unique,
  "Loxodontini" = tree.Atlantogenata.table %>% filter(label=="Loxodontini") %>% pull(node) %>% unique,
  "Loxodona" = tree.Atlantogenata.table %>% filter(label=="Loxodona") %>% pull(node) %>% unique,
  "Elephantidae" = tree.Atlantogenata.table %>% filter(label=="Elephantidae") %>% pull(node) %>% unique,
  "Elephantina" = tree.Atlantogenata.table %>% filter(label=="Elephantina") %>% pull(node) %>% unique,
  "Mammuthus" = tree.Atlantogenata.table %>% filter(label=="Mammuthus") %>% pull(node) %>% unique,
  "Elephantimorpha" = tree.Atlantogenata.table %>% filter(label=="Elephantimorpha") %>% pull(node) %>% unique,
  "Proboscidea" = tree.Atlantogenata.table %>% filter(label=="Proboscidea") %>% pull(node) %>% unique,
  "Tethytheria" = tree.Atlantogenata.table %>% filter(label=="Tethytheria") %>% pull(node) %>% unique,
  "Paenungulata" = tree.Atlantogenata.table %>% filter(label=="Paenungulata") %>% pull(node) %>% unique,
  "Pseudoungulata" = tree.Atlantogenata.table %>% filter(label=="Pseudoungulata") %>% pull(node) %>% unique,
  "Afrotheria" = tree.Atlantogenata.table %>% filter(label=="Afrotheria") %>% pull(node) %>% unique,
  "Afrosoricida" = tree.Atlantogenata.table %>% filter(label=="Afrosoricida") %>% pull(node) %>% unique,
  "Macroscelidae" = tree.Atlantogenata.table %>% filter(label=="Macroscelidae") %>% pull(node) %>% unique,
  "Xenarthra" = tree.Atlantogenata.table %>% filter(label=="Xenarthra") %>% pull(node) %>% unique,
  "Atlantogenata" = tree.Atlantogenata.table %>% filter(label=="Atlantogenata") %>% pull(node) %>% unique
)

atlantogenata.name.node.list <- names(tree.Atlantogenata.nodes) %>% set_names(tree.Atlantogenata.nodes)

tree.Atlantogenata.table %<>% mutate(label = sapply(seq_along(node), function(x, l=atlantogenata.name.node.list, Label=label, Node=as.character(node)){ifelse(Node[x] %in% names(l), l[[Node[x]]], Label[x])}))

species.highlight <- c("Loxodonta africana", "Trichechus manatus", "Procavia capensis", "Orycteropus afer", "Dasypus novemcinctus", "Choloepus hoffmanni", "Chrysochloris asiatica", "Echinops telfairi", "Elephantulus edwardii") %>% 
  set_names(.,.)

tree.Atlantogenata.table %<>% mutate(
  highlight = sapply(seq_along(Species), function(x, l=species.highlight, Label=Species, PhyloPicUID=phyloPicUID){ifelse(Label[x]%in% names(l), PhyloPicUID[x], NA)})#,
  # size.change = node %>% sapply(., branch.size.ratio, tree.data=tree.Atlantogenata.table)
)


tree.Atlantogenata <- tree.Atlantogenata.table %>% as.treedata()


st.tree %>% as_tibble %>% write_csv(path = outTree.f)
tree.Atlantogenata.table %>% write_csv(path=outTree.Atlantogenata.f)
