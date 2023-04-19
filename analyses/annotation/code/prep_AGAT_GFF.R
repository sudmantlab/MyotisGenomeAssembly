library(tidyverse)

#in.gff = '../data/toga/mMyoLuc1/hg38/_query_annotation.agat.gff'
in.gff = snakemake@input[['gff']]
#in.geneids = '../data/toga/mMyoLuc1/hg38/query_isoforms.tsv'
in.geneids = snakemake@input[['geneids']]
#out.gff = '../data/toga/mMyoLuc1/hg38/query_annotation.agat.asMRNA.gff3'
out.gff = snakemake@output[[1]]

write("Reading files", stderr())

gff = read_tsv(in.gff, skip=1, col_names=F)
geneids = read_tsv(in.geneids)

write("Getting ID, Parent, and Name", stderr())
gff.2 <- gff %>% mutate(ID = X9 %>% str_extract("ID=[A-Za-z0-9]+;")%>% str_remove("ID=") %>% str_remove(";"), Parent = X9 %>% str_extract("Parent=[A-Za-z0-9.]+")%>% str_remove("Parent=") %>% str_remove(";"), Name = X9 %>% str_extract("Name=[A-Za-z0-9.]+")%>% str_remove("Name=") %>% str_remove(';')) %>% select(-X9)

write("Merging GFF with gene spans by Name", stderr())
gff.3 <- gff.2 %>% left_join(geneids, by=c('Name'='Projection_ID'))

write("Using the Region ID (gene span name) to fill in missing parents", stderr())
gff.4 <- gff.3 %>% mutate(Parent= sapply(1:length(Parent), function(x, P=Parent, R=Region_ID) { ifelse(is.na(P[x]), R[x], P[x]) }))

# When a line is a gene, it won't have a region ID. 
noregion.IDs <- gff.4[gff.4$X3=='gene',] %>% pull(ID)

write("Genes are really mRNAs - fix that", stderr())
gff.4[!is.na(gff.4$Region_ID),'X3'] <- 'mRNA'

write("Dropping fake UTRs", stderr())
gff.5 <- gff.4 %>% filter(!str_detect(X3, 'UTR'))

write("Keep mRNAs without a parent gene as a parent gene, and give them a unique name. AGAT will add an mRNA line later.", stderr())
gff.5[gff.5$X3=='gene',] <- gff.5[gff.5$X3=='gene',] %>% unite('ID', c('ID','Name'), sep="-", remove=F)

write("Fix parents of these orphan transcripts.", stderr())
gff.5b <- gff.5[gff.5$Parent %in% noregion.IDs,] %>% left_join(gff.4[gff.4$X3=='gene',] %>% unite('Parent2', c('ID','Name'), sep="-", remove=F) %>% select(ID,Parent2), by=c('Parent'='ID')) %>% select(-Parent)

gff.6 <- left_join(gff.5, gff.5b)

gff.7 <- gff.6 %>% mutate(Parent = sapply(1:length(Parent), function(n, P=Parent, Q=Parent2){ifelse(is.na(Q[n]),P[n], Q[n])}))

write("Add back in labels to each element of the attributes column.", stderr())
gff.8 <- gff.7 %>% 
	mutate(
		ID = str_c('ID=', ID), 
		Parent= str_c('Parent=', Parent), 
		Name = sapply(Name, function(x){ifelse(x %>% is.na, yes='', no=paste0('Name=', x))})
		)
		
write("Unite the elements of the attributes column.", stderr())
gff.9 <- gff.8 %>% 
	unite('X9', c('ID','Parent','Name'), sep=';') %>% 
	select(-Region_ID, -Parent2) %>% 
	mutate(X9=X9 %>% str_remove(';$') %>% str_remove_all(';NA')) 

write("Write the file...", stderr())
write(class(out.gff), stderr())
gff.9 %>% 
	write_tsv(out.gff, col_names=F)

write("Done!", stderr())
