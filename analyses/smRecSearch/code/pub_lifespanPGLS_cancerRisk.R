#!/usr/bin/env Rscript

library(tidyverse)
library(ape)
library(treeio)
library(tidytree)
library(ggtree)
library(nlme)
library(broom.mixed)

clade.highlight <- c("Loxodontini", "Loxodontaafricana", "Loxodona", 
                     "Loxodonta cyclotis", "Palaeoloxodon antiquus", 
                     "Elephantidae", "Elephantina", "Elephas maximus",
                     "Mammuthus", "Mammuthus primigenius", "Mammuthus columbi", 
                     "Elephantimorpha",
                     "Proboscidea", "Mammut americanum", 
                     "Tethytheria", "Trichechus manatus", 
                     "Paenungulata", "Procavia capensis", 
                     "Pseudoungulata", "Orycteropus afer", 
                     "Macroscelidae", "Elephantulus edwardii", 
                     "Afrosoricida", "Chrysochloris asiatica", "Echinops telfairi", 
                     "Afrotheria", 
                     "Xenarthra", "Dasypus novemcinctus", "Choloepus hoffmanni", 
                     "Atlantogenata") %>% str_replace_all(" ", ".")

# anage.f <- "../data/AnAge/anage_build14.txt"
# tree.Atlantogenata.f <- "../output/pubFiles/bodysize-Atlantogenata.csv"
# outFile.cancerSuccep.f <- "../output/pubFiles/cancerRisk-Atlantogenata.csv"
# outFile.table.f <- "../output/pubFiles/lifespan-PGLS.tex"
# outFile.htmltable.f <- "../output/pubFiles/lifespan-PGLS.html"

anage.f <- snakemake@input[["anage"]]
tree.Atlantogenata.f <- snakemake@input[["Atlantogenata"]]
outFile.cancerSuccep.f <- snakemake@output[["cancerSuccep"]]
outFile.table.f <- snakemake@output[["latex"]]
outFile.htmltable.f <- snakemake@output[["html"]]


anage_full <- read_tsv(anage.f, col_types = list("Weaning (days)"=col_number(), "Weaning weight (g)"=col_number(), "References" = col_character())) %>% unite("label", Genus, Species, sep=".")
anage <- anage_full %>% 
  select(label,Family, "Maximum longevity (yrs)", "Adult weight (g)") %>% 
  mutate(lnSize.anage = sapply(
    `Adult weight (g)`, 
    function(x){ifelse(is.na(x), NA, log(x))}
  ),
  Lifespan = sapply(
    `Maximum longevity (yrs)`, 
    function(x){ifelse(is.na(x), NA, log(x))}
  )
  ) %>% 
  select(label, Family, lnSize.anage, Lifespan)

dir.create(dirname(outFile.cancerSuccep.f), showWarnings = F, recursive = T)

tree.Atlantogenata.table <- read_csv(tree.Atlantogenata.f) %>% 
  mutate(label=str_replace_all(label, " ", "."))

class(tree.Atlantogenata.table) %<>% str_replace("spec_tbl_df", "tbl_tree")

tree.Atlantogenata <- tree.Atlantogenata.table %>% as.treedata()

tree.size.lifespan.Atlantogenata.table <- left_join(
  # tree.Atlantogenata %>% as_tibble %>% select(parent, node, branch.length, label, lnSize), 
  tree.Atlantogenata.table,
  anage %>% select(-lnSize.anage), 
  by="label"
)

tree.size.lifespan.Atlantogenata <- tree.size.lifespan.Atlantogenata.table %>% 
  as.treedata

#### Make the model

atl.noLife <- tree.size.lifespan.Atlantogenata.table %>%
  filter(is.na(Lifespan)) %>% 
  pull(label)

atl.withLife <- tree.size.lifespan.Atlantogenata.table %>% 
  filter(!is.na(Lifespan)) %>% 
  pull(label)

tree.size.lifespan.Atlantogenata.withLife <- tree.size.lifespan.Atlantogenata %>%
  treeio::drop.tip(atl.noLife) %>% 
  as.phylo()

dat.size.lifespan.Atlantogenata.withLife <- tree.size.lifespan.Atlantogenata.table %>% 
  filter(label %in% atl.withLife) %>% 
  as.data.frame()

rownames(dat.size.lifespan.Atlantogenata.withLife) <- dat.size.lifespan.Atlantogenata.withLife$label

dat.slawl.min <- dat.size.lifespan.Atlantogenata.withLife %>% select(label, Lifespan, lnSize)
rownames(dat.slawl.min) <- dat.slawl.min$label

size.lifespan.Atlantogenata.lm <- gls(
  Lifespan ~ lnSize,
  data=dat.slawl.min,
  correlation = corBrownian(1, tree.size.lifespan.Atlantogenata.withLife, form=~lnSize | label)
)

### Predict new lifespans

d.lifespan.atlantogenata.long <- tree.size.lifespan.Atlantogenata %>% 
  as_tibble() %>% 
  dplyr::select(node, lnSize = Median, lnSize.low = `95%CI_Low`, lnSize.high = `95%CI_High`) %>% 
  pivot_longer(-node, values_to="lnSize")

d.lifespan.atlantogenata.long.predict <-  bind_cols(
  d.lifespan.atlantogenata.long %>% dplyr::select(-lnSize),
  lnSize=predict(
    object = size.lifespan.Atlantogenata.lm, 
    newdata = d.lifespan.atlantogenata.long
  ) %>% exp
) %>% 
  mutate(name = str_replace(name, "lnSize", "Lifespan"))

d.lifespan.atlantogenata <- 
  bind_rows(
    d.lifespan.atlantogenata.long, 
    d.lifespan.atlantogenata.long.predict
  ) %>% 
  pivot_wider(id_cols = node, names_from=name, values_from = lnSize)

tree.size.lifespan.Atlantogenata.full.table <- tree.size.lifespan.Atlantogenata %>% 
  as_tibble %>% 
  select(-Lifespan, -lnSize, -`95%CI_Low`, -`95%CI_High`) %>% 
  left_join(d.lifespan.atlantogenata, by='node') %>% 
### Cancer Susceptibility
  mutate(
    K1 = (Lifespan)^6 * exp(lnSize),
    lnK1 = 6*log(Lifespan) + lnSize
  )

cancerSucceptibility.Atlantogenata.table <- tree.size.lifespan.Atlantogenata.full.table %>% 
  mutate(
    Ancestor = label %>% 
      sapply(
        .,
        function(n, d=tree.size.lifespan.Atlantogenata.full.table){
          if(is.na(n)){return(NA)}
          #print(n)
          #print(length(n))
          p.node <- d %>% filter(label==n) %>% pull(parent) %>% unique
          d %>% filter(node == p.node) %>% pull(label) %>% unique %>% return
        }
      ) %>% str_replace_all("\\.", " ") %>% 
      str_replace("Root", "Atlantogenata"),
    K2 = parent %>% 
      sapply(
        .,
        function(n, d=tree.size.lifespan.Atlantogenata.full.table){
          d %>% filter(node==n) %>% pull(K1) %>% unique
        }
      ),
    lnK2 = parent %>% 
      sapply(
        .,
        function(n, d=tree.size.lifespan.Atlantogenata.full.table){
          d %>% filter(node==n) %>% pull(lnK1) %>% unique
        }
      ),
    CancerSucceptabilityChange = K2/K1,
    log2CancerSucceptabilityChange = log2(K2)-log2(K1),
    label = label %>% str_replace_all("\\.", " ")
  )

### Stats for PGLS:


### Outputs

# print("hi")

cancerSucceptibility.Atlantogenata.table %>% write_csv(path=outFile.cancerSuccep.f)

a <- capture.output(
  stargazer::stargazer(
    size.lifespan.Atlantogenata.lm, type = "latex", 
    title="PGLS: ln(Lifespan) ~ ln(Size)", dep.var.labels = "ln(Lifespan)", 
    covariate.labels = "ln(Size)", colnames = F, dep.var.caption="", 
    out = outFile.table.f, header=F
  )
)
a <- capture.output(
  stargazer::stargazer(
    size.lifespan.Atlantogenata.lm, type = "html", 
    title="PGLS: ln(Lifespan) ~ ln(Size)", dep.var.labels = "ln(Lifespan)", 
    covariate.labels = "ln(Size)", colnames = F, dep.var.caption="", 
    out = outFile.htmltable.f, header=F
  )
)
