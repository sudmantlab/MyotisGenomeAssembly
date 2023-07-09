library(tidyverse)

heart <- read_tsv('output/rsem/Myotis_velifer/mMyoVel1-gff_final_curated-finalAnnotation/LE001-heart.genes.results') %>% 
           mutate(tissue = 'heart', genome='mMyoVel1', anno='gff_final_curated')
lung <- read_tsv('output/rsem/Myotis_velifer/mMyoVel1-gff_final_curated-finalAnnotation/LE001-lung.genes.results') %>%
           mutate(tissue = 'lung', genome='mMyoVel1', anno='gff_final_curated')
brain <- read_tsv('output/rsem/Myotis_velifer/mMyoVel1-gff_final_curated-finalAnnotation/LE001-brain.genes.results') %>%
           mutate(tissue = 'brain', genome='mMyoVel1', anno='gff_final_curated')
kidney <- read_tsv('output/rsem/Myotis_velifer/mMyoVel1-gff_final_curated-finalAnnotation/LE001-kidney.genes.results') %>%
           mutate(tissue = 'kidney', genome='mMyoVel1', anno='gff_final_curated')
testis <- read_tsv('output/rsem/Myotis_velifer/mMyoVel1-gff_final_curated-finalAnnotation/LE001-testis.genes.results') %>%
           mutate(tissue = 'testis', genome='mMyoVel1', anno='gff_final_curated')
pancreas <- read_tsv('output/rsem/Myotis_velifer/mMyoVel1-gff_final_curated-finalAnnotation/LE001-pancreas.genes.results') %>%
           mutate(tissue = 'pancreas', genome='mMyoVel1', anno='gff_final_curated')

all_tissues <- bind_rows(heart, lung, brain, kidney,testis, pancreas) 

all_tissues_wide_tpm <- all_tissues %>% select(genome, anno, tissue, gene_id, transcript_ids=`transcript_id(s)`, TPM) %>% mutate(tissue = tissue %>% str_c(.,'.TPM')) %>% pivot_wider(names_from=tissue, values_from=TPM)

