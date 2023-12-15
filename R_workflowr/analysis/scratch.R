library(tidyverse)
library(ape)
library(ggtree)

k99Dat = read_csv('/mnt/c/users/manue/OneDrive/Desktop/K99GenomeDat.csv') %>% 
  # select(-LQ) %>% 
  rename(BUSCO=`BUSCO Complete`,
         label=Species)

species_cellLine <- read_tsv('../data/batcelllines.tsv') %>% 
  rename(label = Species) %>% 
  mutate(label = label %>% str_replace_all(" ", "_"))

# species = k99Dat %>% pull(label)
species = species_cellLine %>% pull(label)

tree = ape::read.tree('../data/timetree/Run2_real_combined_allGenesAllSpecies/FigTree.nwk')

p = tree %>% ape::keep.tip(species[species %in% tree$tip.label])

plot.base <- revts(p %>% left_join(k99Dat) %>% ggtree(aes(color=LQ), size=2) + 
        geom_tiplab(color='Black', offset = 0.01) + 
        theme_tree2()) + 
  scale_color_viridis_c(option='C') + 
  scale_x_continuous(labels = c(50,12,8,3,0),
                     breaks = c(-0.5,-0.12,-0.08,-0.3,0),
                     limits=c(-0.55, 0.1)) #+
  # theme(legend.position = c(0.05,0.8))


k99Dat.noLQ <- k99Dat %>% select(-LQ)

plot.base.dat <- plot.base %<+% k99Dat.noLQ


plot.base.dat + 
  geom_facet(panel = "BUSCO", data = k99Dat.noLQ, geom = geom_col, 
             mapping=aes(x = BUSCO), orientation='y', color='black', fill='lightblue')+ 
  geom_facet(panel = "BUSCO", data = k99Dat.noLQ, geom = geom_text, 
             mapping=aes(x=BUSCO, label = BUSCO %>% scales::label_percent(scale = 1)(.)), 
             orientation='y', color='black', hjust=1.2) + 
  geom_facet(panel = "Cell Lines", data = k99Dat.noLQ, geom = geom_col, 
             mapping=aes(x = `Cell Lines`), orientation='y', color='black', fill='lightblue')+ 
  geom_facet(panel = "Cell Lines", data = k99Dat.noLQ, geom = geom_text, 
             mapping=aes(x=`Cell Lines`, label = `Cell Lines`), hjust=1.2, orientation='y', color='black') +
  geom_facet(panel = "Paired\nTissue", data = k99Dat.noLQ, geom = geom_col, 
             mapping=aes(x = `Paired Tissue`), orientation='y', color='black', fill='lightblue')+ 
  geom_facet(panel = "Paired\nTissue", data = k99Dat.noLQ, geom = geom_text, 
             mapping=aes(x=`Paired Tissue`, label = `Paired Tissue`), hjust=1.2, orientation='y', color='black') + 
  
  xlim_tree(c(-0.15,0.4))
  
p.celllines <- species_cellLine %>% 
  arrange(desc(Count)) %>% 
  mutate(label = label %>% str_replace("_", ' ')) %>% 
  mutate(label=factor(label, levels=label %>% unique)) %>% 
  ggplot(
    aes(
      x=Count,
      y=label,
      fill=Clade
    )
  ) + 
    geom_bar(stat='identity') + 
    geom_text(aes(label=Count), hjust=1, fontface='bold') + 
  scale_fill_brewer('Clade', palette = 'Dark2') + 
  labs(y='',
       x='# Cell Lines') + 
  theme_pubr() + 
  labs_pubr() + 
  theme(
    axis.text.y = element_text(face = 'bold.italic'),
    legend.position = c(0.75,0.75)
  )

p.celllines  
