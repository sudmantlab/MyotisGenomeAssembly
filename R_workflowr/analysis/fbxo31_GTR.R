tr.fbxo31 <- read.newick("~/projects/fbxo31/output/generax/FBXO31/correction/results/FBXO31/geneTree.newick") %>% as_tibble() %>% separate(label, c('Gene', 'Copy_Species'), sep = '\\.', remove=F) %>% separate(Copy_Species, c("Copy", "Species"), sep="_") %>% unite(label_clean, Gene, Copy, sep=".") %>% mutate(Species = Species %>% str_replace("Myotis", "Myotis_"), highlight = ifelse(str_detect(Species, "Myotis"), Species, NA)) %>% as.treedata()


l.names <- tr.fbxo31 %>% as_tibble() %>% pull(label.y) %>% set_names(tr.fbxo31 %>% as_tibble() %>% pull(label))

tr.fbxo31@phylo$tip.label <- l.names[tr.fbxo31@phylo$tip.label]

plotly::ggplotly(tr.fbxo31 %>% ggtree(aes(text=label)))

nodes.ortho <- c(
  tr.fbxo31 %>% as_tibble %>% filter(str_detect(Species, "Myotis", negate=T)) %>% pull(node))
  85,84,82,83,81,77,78,79,80
  )


tr.fbxo31 %>% as_tibble() %>% filter(node %in% nodes.ortho)

tr.fbxo31 %>% 
  ggtree(layout = 'circular') + 
  # geom_highlight(mapping=aes(subset=node %in% nodes.ortho), fill='blue') + 
  geom_strip(
    "FBXO31.1_Loxodontaafricana", 'FBXO31.1_Rhinolophusferrumequinum', 
    barsize=1, color='black', offset=1,
    label="FBXO31\n1:1\northologs", offset.text=5
  ) + 
  geom_strip(
    "FBXO31.24_Myotisvolans", 'FBXO31.30_Myotisevotis', 
    barsize=1, color='blue', offset=1,
    label="Myotis\nCanon\nFBXO31", offset.text=5
  ) + 
  geom_strip(
    "FBXO31.33_Myotisoccultus", 'FBXO31.5_Myotislucifugus', 
    barsize=1, color='red', offset = 1,
    label="Recurrent\nMyotis\nFBXO31\nexpansion", offset.text=12
  ) + 
  geom_tippoint(mapping=aes(subset=!is.na(highlight), color=highlight)) +
  scale_color_manual(values=species_color, guide=guide_none()) + 
  xlim(c(0,50))
