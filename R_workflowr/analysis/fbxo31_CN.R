tr.species <- read.tree('../data/stableTraits/UphamEtAl2019_newChiroptera.tree')


dat <- read_tsv('../../../fbxo31/output/myotis_FBPCount.tsv')

dat <- dat %>% filter(!Species %in% setdiff(dat$Species, tr.species$tip.label))

dat.fbxo31 <- dat %>% 
  filter(label == "FBXO31") %>% 
  dplyr::rename(Gene = label) %>% 
  dplyr::rename(
         label = Species)

tr.species <- keep.tip(tr.species, dat.fbxo31$label)

tr.species$tip.label <- str_replace(tr.species$tip.label, "_", " ")

dat.fbxo31$label <- str_replace(dat.fbxo31$label, "_", " ")

dat.fbxo31.anno <- read_tsv('../data/genes/gene_copy_MyotisAnnotation.tsv') %>% 
  dplyr::filter(human_gene_name == "FBXO31") %>% 
  dplyr::select(-human_gene_name) %>% 
  dplyr::rename(label=species) %>% 
  dplyr::mutate(label = str_replace(label, "_", " "))

our_genomes = c(
  "Myotis_auriculus",
  "Myotis_californicus",
  "Myotis_occultus",
  "Myotis_lucifugus",
  "Myotis_yumanensis",
  "Myotis_volans",
  "Myotis_velifer",
  "Myotis_evotis",
  "Myotis_thysanodes"
)

species_color = ggsci::pal_d3(palette = "category20")(length(our_genomes)+1) %>% 
  set_names(., c(our_genomes, "Other"))
species_color["Myotis_evotis"] = "#17BECFFF"
species_color["Other"] = "#7F7F7FFF"

names(species_color) <- str_replace(names(species_color), "_"," ")

p.tr <- tr.species %>% 
  ggtree() + 
  geom_tiplab(aes(color=label), fontface="bold.italic") + 
  coord_cartesian(clip="off") +
  scale_color_manual(values=c(species_color, "black"), guide=guide_none()) + 
  xlim(c(NA, 75))

p.tr

p.dat <- dat.fbxo31 %>% 
  ggplot(
    aes(
      y=label,
      x=n
    )
  ) + 
  labs(
    x="FBXO31 Copy Number"
  ) +
  geom_segment(
    aes(x=0,xend=n,
        yend=label,
        color = "RBHB")
  ) + 
  geom_point(
    aes(color = "RBHB")
  ) + 
  geom_segment(
    data = dat.fbxo31.anno,
    aes(x=0,xend=n,
        yend=label,
        color = "Full Copy")
  ) + 
  geom_point(
    data = dat.fbxo31.anno,
    aes(color = "Full Copy")
  ) + 
  scale_color_manual("Type", values=c("RBHB"='black', "Full Copy"='goldenrod')) +
  theme_pubr()+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'right'
  )

p.dat

# p.tr %>% aplot::insert_right(p.dat)

p.dat %>% aplot::insert_left(p.tr)

