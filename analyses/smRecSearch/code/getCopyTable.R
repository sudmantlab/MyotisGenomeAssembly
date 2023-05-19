suppressMessages(suppressWarnings(require(dplyr, quietly = T, warn.conflicts = F)))
suppressMessages(suppressWarnings(require(readr, quietly = T, warn.conflicts = F)))
options(readr.num_columns = 0)
suppressMessages(suppressWarnings(require(plyranges, quietly = T, warn.conflicts = F)))
suppressMessages(suppressWarnings(source("code/generalFunctions.R")))

##############################################################
imagesavename <- stringr::str_c(snakemake@output[["copyTable"]],".RData")
save.image(file=imagesavename)

print("Getting Core Count....")
core_count <- snakemake@threads[1]
print("Core Count: ", str(core_count))

print("Reading Translation Table...")
translation.table <- read_tsv(snakemake@input[["TranslationTable"]])
print("Done.")

print("Reading RBB file...")
rbb <- get.rbb(snakemake@input[["RBB_file"]])%>%
  plyranges::mutate(Gene = sapply(QName, function(id, tt = translation.table){tt$Name[which(tt$ID == id)][1]}))
print("Done.")

print("Reading ECNC...")
ecnc <- get.ecnc(snakemake@input[["ECNC_file"]]) %>%
  full_join(., translation.table, by=c("name"="ID"))
print("Done.")

print("Updating save.image with RBB, ECNC...")
save.image(file=imagesavename)
print("Done.")

print("Starting reduce.bed...")
# reduced <- rbb %>%
#   reduce.bed(., core_count=core_count)
reduced <- rbb %>%
  reduce.bed.dt(., core_count=core_count)
print("Done.")

print("Writing exonTable...")
reduced %>% 
  write_tsv(snakemake@output[["exonTable"]])

print("Updating save.image with reduced...")
save.image(file=imagesavename)
print("Done.")

print("Starting get.gene.transcript.loci...")
# loci <- reduced %>%
#   get.gene.transcript.loci(., core_count=core_count)
loci <- reduced %>%
  get.gene.transcript.loci.dt(., core_count=core_count)
print("Done.")

print("Updating save.image with loci...")
save.image(file=imagesavename)
print("Done.")

print("Generating finalized BED file")
# final.bed <- get.gene.loci(loci, reduced, core_count=core_count)
# final.bed <- get.gene.loci.dt(loci, reduced, core_count=core_count)
final.bed <- get.gene.loci.dt2(loci, reduced, core_count=core_count)
print("Done.")

print("Updating save.image with final BED...")
save.image(file=imagesavename)
print("Done.")

print("Creating loci.table, loci.max, & loci.ecnc")

loci.table <- loci %>% 
  unite("GeneLocus", Gene, locus, sep="_", remove = F) %>% 
  mutate(TranscriptID = sapply(TranscriptID, . %>% paste(., collapse=",")))

loci.max <- loci %>%
  select(-TranscriptID) %>%
  group_by(Gene) %>%
  summarize(Copy=max(locus)) %>%
  ungroup

loci.ecnc <- ecnc %>%
  select(Name, ECNC) %>%
  group_by(Name) %>%
  summarize(ECNC=max(ECNC)) %>%
  ungroup %>%
  dplyr::rename(Gene = Name) %>%
  full_join(., loci.max, by="Gene")

print("Done.")

print("Writing copyTable")
loci.ecnc %>%
  write_tsv(snakemake@output[["copyTable"]])

loci.table %>%
  as_tibble() %>% 
  write_tsv(snakemake@output[["lociTable"]])

final.bed %>% 
  # write_bed(snakemake@output[["finalBed"]])
  write_tsv(snakemake@output[["finalBed"]])

loci.ecnc %>%
  filter(Copy>1, ECNC>=1.5) %>%
  pull(Gene) %>%
  unique %>%
  write_lines(snakemake@output[["duplicateGeneList"]])

print('Done with everything!')