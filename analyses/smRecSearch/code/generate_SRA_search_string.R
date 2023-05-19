#!/usr/bin/env Rscript

require(tidyverse, quietly = T)

port.table <- read_csv("data/portTable.csv", comment = '"#"', col_names = c("TranslatedPort", "UntranslatedPort", "Genome", "CommonName", "Species"))

special.species <- c(
"Mus musculus",
"Homo sapiens",
"Danio rerio",
"Drosophila melanogaster",
"Saccharomyces cerevisiae",
"Rattus norvegicus",
"Bos taurus",
"Sus scrofa",
"Gallus gallus",
"Macaca mulatta",
"Pan troglodytes",
"Ovis aries",
"Canis lupus",
"Arabidopsis thaliana",
"Papio anubis",
"Macaca fascicularis",
"Xenopus tropicalis",
"Meleagris gallopavo",
"Oryctolagus cuniculus",
"Solanum lycopersicum"
)


species <- port.table %>% pull(Species) %>% unique %>% setdiff(., special.species) %>% str_c(rep('"'), ., rep('"[Organism]')) %>% paste(collapse = " OR ") %>% paste0("(", ., ")")

special.species.str <- special.species %>% str_c(rep('"'), ., rep('"[Organism]')) %>% paste(collapse = " OR ") %>% paste0("(", ., ")")

species.not.special <- str_c(species, " NOT ", special.species.str)

platform <- '"illumina"[Platform]'

access <- '"public"[Access]'

strategy <- '("est"[Strategy] OR "rna seq"[Strategy])'

cat("\n\nMOST SPECIES:\n")

sprintf("((((%s) AND (%s)) AND (%s)) AND (%s))\n", species.not.special, platform, access, strategy) %>% cat

cat("\n\n SPECIAL SPECIES:\n")

for (spc in special.species){
  sprintf('(((("%s"[Organism]) AND (%s)) AND (%s)) AND (%s))\n', spc, platform, access, strategy) %>% cat
}




