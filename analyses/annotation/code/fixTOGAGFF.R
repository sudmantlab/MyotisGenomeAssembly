library(tidyverse)

fixTOGAGFF <- function(gff_in){
	a <- read_tsv(gff_in, col_names=F, skip=1)
	exons <- a %>% filter(X3 != "gene") %>% mutate(ID = X9 %>% str_extract("ID=[A-Za-z0-9]+;")%>% str_remove("ID=") %>% str_remove(";"), Parent = X9 %>% str_extract("Parent=[A-Za-z0-9.]+")%>% str_remove("Parent=") %>% str_remove(";")) %>% select(-X9)
	genes <- a %>% filter(X3 == "gene") %>% mutate(ID = X9 %>% str_extract("ID=[A-Za-z0-9]+;")%>% str_remove("ID=") %>% str_remove(";"), Name = X9 %>% str_extract("Name=[A-Za-z0-9.]+;")%>% str_remove("Name=") %>% str_remove(";"))
	genes.tt <- select(genes, ID, Name)
	exons.renamed <- exons %>% full_join(genes.tt, by=c("Parent" = "ID")) %>% 
		mutate(Parent = Name, ID = str_c(ID, ".", Name)) %>% select(-Name) %>%
		mutate(ID = str_c("ID=",ID), Parent = str_c("Parent=",Parent)) %>% 
		unite("X9", c(ID, Parent), sep=";")
	genes.renamed = genes %>% select(-ID) %>% mutate(X9 = str_c("ID=", Name))
	mrnas <- genes.renamed %>% mutate(X3='mRNA', X10 = X9 %>% 
		str_replace('ID=', 'Parent=TOGA.')) %>% unite("X9", c(X9, X10), sep=';') %>% 
		select(-Name)
	genes.renamed <- genes %>% select(-ID) %>% mutate(X9 = str_c("ID=TOGA.", Name)) %>% select(-Name)
	gff <- bind_rows(genes.renamed, mrnas) %>% bind_rows(., exons.renamed) %>% mutate(X3 = factor(X3, levels = c("gene", "mRNA", "exon", "CDS"))) %>% arrange(X1,X3,X9,X4)
}	

fixTOGAGFF(snakemake@input[[1]]) %>% write_tsv(snakemake@output[[1]], col_names=F)

