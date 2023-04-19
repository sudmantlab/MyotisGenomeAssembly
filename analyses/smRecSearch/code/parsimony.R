suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(ape)))
suppressMessages(suppressWarnings(require(treeio)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(magrittr)))
source("code/generalFunctions.R")
setDTthreads(snakemake@wildcards[["threads"]])


sep <- function(..., suffixes=c(".min", ".max")) {
  # https://stackoverflow.com/questions/42464250/apply-tidyrseparate-over-multiple-columns
  dots <- list(...)
  #n <- stringr::str_count(dots[[1]][[dots[[2]]]], "\\d+")
  separate_(..., into = sprintf("%s_col%d", dots[[2]], suffixes), sep = ",", fill="right")
}

MPR_advanced <- function(name, trait, phylo.unrooted, outgroup){
  dt = MPR(x=trait, phy = phylo.unrooted, outgroup) %>% 
    as.data.table(keep.rownames = T)
  g.min <- paste(name, "min", sep = ".")
  g.max <- paste(name, "max", sep = ".")
  setnames(dt, "lower", g.min)
  setnames(dt, "upper", g.max)
  setnames(dt, "rn", "node")
  # dt[, node:=as.numeric(node)]
  # dt[, (g.min):=as.double(g.min)]
  # dt[, (g.max):=as.double(g.max)]
  setkey(dt, "node")
  melt(dt, id.vars="node")
  }

tree <- ape::read.tree(snakemake@input[["tree"]])
#tree <- read.tree("data/stableTraits/eutheria.tree")


tree.nodenames <- tree$node.label %>% set_names(seq_along(.))
tree$node.label <- tree.nodenames %>% names
# afrotheria.tree <- tree %>% 
#   tree_subset(node=tidytree::MRCA(tree, "Loxodonta.africana", "Dasypus.novemcinctus"), levels_back=0)


geneCopyTable.wide <- fread(snakemake@input[["wideTable"]])
#geneCopyTable.wide <- fread("output/geneCopyTable/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/consolidated-GeneCopyTable_RBB_wide.csv")[label %in% tree$tip.label]
#afrotheria.tree <- drop.tip(afrotheria.tree, tip = tree$tip.label[!tree$tip.label %in% geneCopyTable.wide$label])
# afrotheria.tree.data <- afrotheria.tree %>% 

save.image("debug_parsimony.RData")

tree.data <- tree %>% 
  full_join(geneCopyTable.wide, by="label")


t <- unroot(tree)

g <- tree.data@data %>% select(-node)
#MPR(x=g[,1][[1]], phy =t, "Dasypus.novemcinctus") %>% apply(., 1, paste, collapse=",")
o <- lapply(names(g), function(n, gdt=g, p=t, og="Dasypus.novemcinctus"){MPR_advanced(name=n, trait=gdt[[n]], phylo.unrooted = t, outgroup = og)}) %>% 
  #lapply(melt, id.vars="node") %>% 
  rbindlist() %>% 
  dcast(node~variable, value.var="value") # all this cause cbindlist hasn't been written yet


geneCopyTable.wide <- merge(geneCopyTable.wide, geneCopyTable.wide, by="label", suffixes = c(".min", ".max")) # a hack if I ever saw one

tree.data <- tree %>%
  full_join(geneCopyTable.wide, by="label") %>% 
  as_tibble

save.image("debug_parsimony.RData")

o[, node:=as(node, tree.data$label %>% class)]

tree.table <- coalesce_join(tree.data, o, by=c("label"="node"))

tree.parsimonious <-tree.table %>% as.treedata()

tree.parsimonious@phylo$node.label %<>% sapply(function(x, l=tree.nodenames){l[[x]]}, USE.NAMES=F)

write.beast(tree.parsimonious, snakemake@output[[1]])



