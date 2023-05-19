suppressMessages(suppressWarnings(require(dplyr, quietly = T, warn.conflicts = F)))
suppressMessages(suppressWarnings(require(readr, quietly = T, warn.conflicts = F)))

suppressMessages(suppressWarnings(source("code/generalFunctions.R")))

warning("DEPRECIATED SCRIPT!!!!")
##############################################################
imagesavename <- stringr::str_c(basename(snakemake@output[["duplicateGeneList"]]),".RData")
save.image(file=imagesavename)

print("Getting Core Count....")
core_count <- snakemake@threads[1]
print("Core Count: ", str(core_count))

print("Reading Translation Table...")
translation.table <- read_tsv(snakemake@input[["TranslationTable"]]) %>%
  filter_func
print("Done.")

print("Reading RBB file...")
rbb <- get.rbb(snakemake@input[["RBB_file"]])%>%
  mutate(Gene = sapply(QName, function(id, tt = translation.table){tt$Name[which(tt$ID == id)][1]})) %>%
  filter(Gene %in% translation.table$Name)
print("Done.")

print("Reading ECNC...")
ecnc <- get.ecnc(snakemake@input[["ECNC_file"]]) %>%
  right_join(., translation.table, by=c("name"="ID"))
print("Done.")

print("Updating save.image with RBB, ECNC...")
save.image(file=imagesavename)
print("Done.")

print("Starting reduce.bed...")
reduced <- rbb %>%
  reduce.bed(., core_count=core_count)
print("Done.")

print("Updating save.image with reduced...")
save.image(file=imagesavename)
print("Done.")

print("Starting get.gene.transcript.loci...")
loci <- reduced %>%
  get.gene.transcript.loci(., core_count=core_count)
print("Done.")

print("Updating save.image with loci...")
save.image(file=imagesavename)
print("Done.")

print("Creating loci.max & loci.ecnc")
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

print("Outputting final files...")

reduced %>% 
  as_tibble() %>% 
  write_tsv(snakemake@output[["exonTable"]])

loci.ecnc %>%
  write_tsv(snakemake@output[["copyTable"]])

loci %>%
  as_tibble() %>% 
  write_tsv(snakemake@output[["lociTable"]])

loci.ecnc %>%
  filter(Copy>1, ECNC>=1.5) %>%
  pull(Gene) %>%
  unique %>%
  write_lines(snakemake@output[["duplicateGeneList"]])

print('Done with everything!')
