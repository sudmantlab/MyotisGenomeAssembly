suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
options(readr.num_columns = 0)
source("code/generalFunctions.R")

genomeTable <- read_csv(snakemake@input[["portTable"]], comment = '"#"', col_names = c("TPort", "UPort", "Genome", "Common Name", "Species")) %>% 
  filter(!is.na(Species), Genome != "proCap3") %>% 
  select(-TPort, -UPort) %>% 
  mutate(
    "Species.strict" = str_split(Species, " ") %>% sapply(., function(x){paste(x[1], x[2], sep=" ")}), 
    "label" = str_split(Species.strict, " ") %>% sapply(., function(x){paste(x[1], x[2], sep=".")}),
    "Common Name" = str_remove(`Common Name`, "\\([A-Za-z]+\\)|West Indian"),
    "Spc.short" = str_split(Species, " ") %>% sapply(., function(x){paste(x[1] %>% str_trunc(2, "right", "."), x[2], sep=" ")})
  )%>%
  group_by(label) %>%
  summarize(
    Genome = toString(Genome),
    Species.strict = first(Species.strict),
    `Common Name` = first(`Common Name`),
    Species = first(Species),
    Spc.short = first(Spc.short),
    BestGenome=str_split(Genome, ", ") %>% sapply(. %>% tail(n=1))
  ) %>%
  ungroup() %>% 
  distinct() %>% 
  arrange(label) %>% 
  mutate(phyloPicUID = phylopic_uid_vector(str_replace_all(Species, " ", "_")))


genomeTable[genomeTable$label=="Orycteropus.afer",]$Genome <- "oryAfe1, oryAfe2"
genomeTable[genomeTable$label=="Orycteropus.afer",]$BestGenome <- "oryAfe2"

genomeTable %>% write_csv(snakemake@output[["genomeTable"]])