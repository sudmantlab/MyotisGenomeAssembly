###############           QC Functions
qc.name.bed <- function(bed){
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  options(readr.num_columns = 0)
  bed %>% 
    separate(name, c("name", "hit"), sep="_") %>% 
    group_by(name) %>% 
    mutate(hit=row_number(hit)-1) %>% 
    ungroup %>% 
    unite(name, hit, col="name")
}


filter_func <- function(df){
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  options(readr.num_columns = 0)
  df %>%
    filter(
      !is.na(Name),  # No unnamed Names, these are Namerally uncharacterized "Names"
      !str_detect(Name, "C[0-9XY]{1,2}orf\\d+"),  # No C-Orfs
      !str_detect(Name, "^LOC\\d+"),   # No uncharacterized loci
      !str_detect(Name, "^LINC\\d+"),   # No Long "Non-Coding" RNA
      !str_detect(Name, "^HLA|^HCP5|^CD74"), # No HLAs
      !str_detect(header, "[Pp]rotocadherin"), # No protocadherins
      !str_detect(Name, "^KRT\\d"),            # No Keratins
      !str_detect(Name, "^H[1-4]|^HIST\\d"), # No Histones
      !str_detect(Name, "^ZNF\\d+"),   # No Zinc Fingers
      !str_detect(Name, "^RP[LS]\\d+"),   # No Riboproteins
      !str_detect(Name, "^OR\\d+"),   # No Olfactory Receptors
      !str_detect(header, "[uU]ncharacterized"),  # No Uncharacterized proteins
      !str_detect(header, "[pP]utative"),  # No putative proteins
      !str_detect(header, "[pP]olyprotein|ERVK|HERV"),  # No viral proteins
      !str_detect(header, "[Ff]ragment")  # No fragment proteins
    )
}


multi_join <- function(list_of_loaded_data, join_func, ...){
  
  suppressMessages(suppressWarnings(require(dplyr, quietly=T)))
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  
  return(output)
}


##############         READ FUNCTIONS

get.rbb <- function(rbb.file){
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  options(readr.num_columns = 0)
  suppressMessages(suppressWarnings(require(rtracklayer, quietly=T)))
  extracolnames <- c("QName" = "character", "QStart" = "character", "QEnd" = "character", "QCov" = "character")
  import.bed(rbb.file, extraCols=extracolnames)
}


get.ecnc <- function(ecnc.file){
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  options(readr.num_columns = 0)
  read_tsv(ecnc.file, col_names = c("name", "ECNC"))
}



############        WORKER FUNCTIONS


bed_blocks_to_rows <- function(bed.df){
  suppressMessages(suppressWarnings(require(magrittr, quietly=T)))
  suppressMessages(suppressWarnings(require(plyranges, quietly=T)))
  bed.df %>%
    mutate(block_count=sapply(blocks, length)) %>%
    expand_ranges(blocks, .id="Exon") %>%
    narrow(., start=start(.$blocks), end=end(.$blocks)) %>%
    group_by(Exon) %>%
    mutate(Exon=row_number(Exon)) %>%
    ungroup()
}

bed_rows_to_bed_block <- function(bed.oneGene){
  suppressMessages(suppressWarnings(require(GenomicRanges, quietly = T, warn.conflicts = F)))
  newStrand <- bed.oneGene %>% strand %>% unique
  new.GRange <- GRanges(
    seqnames = bed.oneGene %>% seqnames %>% unique,
    ranges = IRanges(
      start = bed.oneGene %>% start %>% min,
      end = bed.oneGene %>% end %>% max
    ),
    strand = if(length(newStrand==1)) newStrand else Rle(strand("*")), 
    name = bed.oneGene$name %>% unique, 
    score = 1, 
    itemRgb = "#FF0000",
    thick = IRanges(
      start = bed.oneGene %>% start %>% min,
      end = bed.oneGene %>% end %>% max
    ),
    blocks = bed.oneGene %>% ranges %>% shift(-(bed.oneGene %>% start %>% min - 1)) %>% IRangesList
  )
}


reduce_list_sets <- function(list_lists){
  unique(
    lapply(
      list_lists,
      function(x) {
        unique(
          unlist(
            list_lists[sapply(list_lists, function(y){any(x %in% y)})]
          )
        )
      }
    )
  )
}


reduce.chr.bed <- function(bedrec){
  #bedrec %T>% debugprint1("Starting bed_blocks_to_rows") %>%
  # bed_blocks_to_rows %T>% debugprint1("Starting reduce_ranges_directed") %>% 
  bedrec %>% 
    bed_blocks_to_rows %>% 
    reduce_ranges_directed(
      TranscriptID=unique(name), 
      Gene=unique(Gene), 
      Exon=unique(Exon)
    )
}


get.gene.locus <- function(x, bed) {
  bed %>% 
    filter(TranscriptIDPrime %in% x$TranscriptID[[1]]) %>% 
    mutate(name=x$Gene) %>% 
    bed_rows_to_bed_block
}


collapse.one.locus <- function(df){
  suppressMessages(suppressWarnings(require(magrittr, quietly=T)))
  suppressMessages(suppressWarnings(require(dplyr, quietly=T)))
  df %>% 
    # pull(TranscriptID) %T>% debugprint2(paste0(.[[1]], " First round reduce\n")) %>% 
    # reduce_list_sets %T>% debugprint2(paste0(.[[1]], " Second round reduce\n")) %>% 
    # reduce_list_sets %T>% debugprint2(paste0(.[[1]], " Done!\n"))
    pull(TranscriptID) %>% 
    reduce_list_sets %>% 
    reduce_list_sets 
}


####################                    MASTER FUNCTIONS


reduce.bed <- function(bed.df, core_count){
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  options(readr.num_columns = 0)
  suppressMessages(suppressWarnings(require(parallel, quietly = T)))
  suppressMessages(suppressWarnings(require(plyranges, quietly=T)))
  print("Starting ForkCluster for parLapplys")
  cl <- parallel::makeForkCluster(core_count, outfile = "")
  print("Exporting needed functions to cluster")
  parallel::clusterExport(cl, c("reduce.chr.bed"))
  print("Reducing bed file")
  reduced.bed <- bed.df %>%
    # Step 1: Expand BED12 into BED6
    # Breaks apart hits into exons in preparation for merging
    #bed_blocks_to_rows %>%
    # Step 2: Combine overlapping exons, while keeping info of the TranscriptID and Gene of origin.
    # Keeps info on exon number just in case
    split(., seqnames(.)) %>%
    clusterApply(
      cl=cl, 
      x = ., 
      fun = reduce.chr.bed
    ) %>%
    Reduce(append, .) %>%
    # Step 3: Remove Many-Genes-to-One Hits at the Exon levels
    # The whole point is that we expect Many-to-One at the Transcript level; what we don't want is exons that have different genes matching them
    # If the hit is real, we'll likely have another exon of the hit that is unique, and so overall the transcript will still be represented
    filter(sapply(Gene, length)<=1) %>%
    # Step 3B: convert the CharacterList column to a Character column for subsetting
    mutate(
      Gene = sapply(Gene, function(x)x[[1]]),
      TranscriptIDPrime = sapply(TranscriptID, "[[", 1)
    )
  print("Stopping Cluster for parLapply")
  stopCluster(cl)
  return(reduced.bed)
}


get.gene.transcript.loci <- function(reduced.bed, core_count, parallel=T){
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  options(readr.num_columns = 0)
  suppressMessages(suppressWarnings(require(parallel, quietly = T)))
  suppressMessages(suppressWarnings(require(plyranges, quietly=T)))
  if(parallel){
    library(magrittr)
    print("Starting ForkCluster for parLapplys")
    cl <- parallel::makeForkCluster(core_count, outfile = "")
    print("Exporting needed functions to cluster")
    parallel::clusterExport(cl, c("collapse.one.locus"))
    print("Creating Gene-Transcript List")
    gene.transcript.list <- reduced.bed %T>% debugprint1("this is a test") %>%
      # Step 1: Get df of Genes with Transcript IDs mergers
      as.data.table %>%  # Without as.data.table, memory usage here is prohibitive
      select(TranscriptID, Gene, TranscriptIDPrime, seqnames) %T>% debugprint1("before split") %>% 
      # Step 2: Get list of transcripts represented, regardless of exon count
      # Lots of redundancies here because we're still on a per-exon level
      # Step 3: Reduce transcript lists into a list of transcripts per Gene Loci
      #split(.,list(.$seqnames, .$Gene), drop = T)
      split(., .$Gene, drop = T) %T>% debugprint2("after split") %>% 
      # print("Merging overlapping transcripts")
      #final.df <- gene.transcript.list %>%
      #this is slow so it needs to be parallelized
      clusterApply(
        cl = cl, 
        x = ., 
        fun = collapse.one.locus
      ) %>%
      enframe %>%
      unnest %>%
      dplyr::rename(Gene="name", TranscriptID="value") %>%
      group_by(Gene) %>%
      mutate(locus = row_number(Gene)) %>% 
      ungroup
    stopCluster(cl)
    print("Stopping Cluster for parLapply")
  } else {
    library(magrittr)
    print("Creating Gene-Transcript List")
    gene.transcript.list <- reduced.bed %>%
      # Step 1: Get df of Genes with Transcript IDs mergers
      #mcols %>%
      #as_tibble %>% 
      as.data.table %>%  # Without as.data.table, memory usage here is prohibitive
      #arrange(Gene, TranscriptID) %>%
      # No longer need exons here
      select(TranscriptID, Gene, TranscriptIDPrime) %>%
      # Step 2: Get list of transcripts represented, regardless of exon count
      # Lots of redundancies here because we're still on a per-exon level
      # distinct %>%
      # Step 3: Reduce transcript lists into a list of transcripts per Gene Loci
      #split(.,list(.$seqnames, .$Gene), drop = T)
      split(., .$Gene, drop = T)%>%
      # print("Merging overlapping transcripts")
      #final.df <- gene.transcript.list %>%
      #this is slow so it needs to be parallelized
      lapply(. %>% 
                  pull(TranscriptID) %T>% print(.[1]) %>%
                  #str_split(., ",") %>% 
                  reduce_list_sets %T>% 
                  reduce_list_sets # %>% 
                #sapply(., . %>% paste(., collapse=","))
      ) %>%
      enframe %>%
      unnest %>%
      dplyr::rename(Gene="name", TranscriptID="value") %>%
      #separate("Gene", c("seqnames", "Gene")) %>% 
      #distinct %>% 
      group_by(Gene) %>%
      mutate(locus = row_number(Gene)) %>% 
      ungroup
  }
  #return(final.df)
  return(gene.transcript.list)
}


get.gene.loci <- function(geneLociTable, bed.rbb, core_count, parallel=T){
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  options(readr.num_columns = 0)
  suppressMessages(suppressWarnings(require(parallel, quietly = T)))
  suppressMessages(suppressWarnings(require(plyranges, quietly=T)))
  print("Creating Gene-Transcript List")
  gene.transcript.list <- geneLociTable %>% 
    unite("Gene", Gene, locus, sep="_") %>% 
    split(., .$Gene)
  if(parallel){
    print("Starting ForkCluster for parLapplys")
    cl <- parallel::makeForkCluster(core_count, outfile = "")
    print("Exporting needed functions to cluster")
    parallel::clusterExport(
      cl,
      c("%>%", "bed_rows_to_bed_block", "mutate", "filter", "%in%")
    )
    print("Merging BED Transcripts")
    final.bed <- gene.transcript.list %>% 
        parLapply(
        cl,
        ., 
        get.gene.locus,
        bed=bed.rbb
      ) %>% 
      bind_ranges
    print("Stopping Cluster for parLapply")
    stopCluster(cl)
  } else {
    library(magrittr)
    print("Merging BED Transcripts")
    final.bed <- gene.transcript.list %>% 
      lapply(
        ., 
        function(x, bed=bed.df) {
          bed %>%
            filter(TranscriptIDPrime %in% x$TranscriptID[[1]]) %>%
            mutate(name=x$Gene) %T>% debugprint(., bla=F) %>%
            bed_rows_to_bed_block %T>% debugprint(., bla=T)
        }
      ) %>% 
      bind_ranges
  }
  return(final.bed)
}



#############################              DT-based Master Functions



reduce.chr.bed.dt <- function(bedrec){
  suppressMessages(suppressWarnings(require(data.table, quietly=T)))
  suppressMessages(suppressWarnings(require(plyranges, quietly=T)))
  #bedrec %T>% debugprint1("Starting bed_blocks_to_rows") %>%
  # bed_blocks_to_rows %T>% debugprint1("Starting reduce_ranges_directed") %>% 
  bedrec %>% 
    bed_blocks_to_rows %>% 
    reduce_ranges_directed(
      TranscriptID=unique(name), 
      Gene=unique(Gene), 
      Exon=unique(Exon),
      itemRgb= "255,0,0" #col2rgb(itemRgb) %>% rowMeans %>% paste(sep=",")
    ) %>% 
    as.data.table %>% 
    .[sapply(Gene, length) == 1] # Step 3: Remove Many-Genes-to-One Hits at the Exon levels
}


reduce.bed.dt <- function(bed.df, core_count){
  suppressMessages(suppressWarnings(require(data.table, quietly=T)))
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  suppressMessages(suppressWarnings(require(parallel, quietly = T)))
  suppressMessages(suppressWarnings(require(plyranges, quietly=T)))
  setDTthreads(core_count)
  print("Starting ForkCluster for parLapplys")
  cl <- parallel::makeForkCluster(core_count, outfile = "")
  print("Exporting needed functions to cluster")
  parallel::clusterExport(cl, c("reduce.chr.bed"))
  print("Reducing bed file")
  reduced.dt <- bed.df %>%
    # Step 1: Expand BED12 into BED6
    # Breaks apart hits into exons in preparation for merging
    #bed_blocks_to_rows %>%
    # Step 2: Combine overlapping exons, while keeping info of the TranscriptID and Gene of origin.
    # Keeps info on exon number just in case
    split(., seqnames(.)) %>%
    clusterApply(
      cl=cl, 
      x = ., 
      fun = reduce.chr.bed.dt
    ) %>% 
    rbindlist()
  # The whole point is that we expect Many-to-One at the Transcript level; what we don't want is exons that have different genes matching them
  # If the hit is real, we'll likely have another exon of the hit that is unique, and so overall the transcript will still be represented
  reduced.dt[, `:=`(Gene = sapply(Gene, function(x)x[[1]]), TranscriptIDPrime = sapply(TranscriptID, "[[", 1))]
  print("Stopping Cluster for parLapply")
  stopCluster(cl)
  return(reduced.dt)
}


get.gene.transcript.loci.dt <- function(reduced.dt, core_count){
  setDTthreads(core_count)
  suppressMessages(suppressWarnings(require(data.table, quietly=T)))
  reduced.dt[,list(TranscriptID=reduce_list_sets(reduce_list_sets(TranscriptID))), by=Gene][, .(locus=seq_len(.N), TranscriptID), by=Gene]
}


get.gene.locus.dt2 <- function(gene.transcripts, bed.dt, core_count) {
  setDTthreads(core_count)
  bed.dt[TranscriptIDPrime %in% gene.transcripts$TranscriptID[[1]], 
         ][, 
           list(
             seqnames = unique(seqnames),
             start = min(start),
             end = max(end),
             strand = unique(strand) %>% if(length(.==1)) . else "*", 
             name = unique(Gene),
             score = 1, 
             itemRgb = itemRgb[1],
             thickStart = min(start),
             thickEnd = max(end),
             blockCount = .N,
             blockStart = paste(start - min(start)+1, collapse=","),
             blockWidth = paste(end - start, collapse=",")
           )
           ]
}


get.gene.loci.dt2 <- function(geneLociTable, reduced.dt, core_count){
  setDTthreads(core_count)
  suppressMessages(suppressWarnings(require(tidyverse, quietly=T)))
  suppressMessages(suppressWarnings(require(data.table, quietly=T)))
  suppressMessages(suppressWarnings(require(plyranges, quietly=T)))
  print("Merging Gene Loci")
  gene.transcript.list <- geneLociTable %>% 
    unite("Gene", Gene, locus, sep="_") %>% 
    split(., .$Gene)
  final.bed <- gene.transcript.list %>% 
    lapply(
      ., 
      get.gene.locus.dt2,
      bed=reduced.dt,
      core_count=core_count
    ) %>% 
    rbindlist
  final.bed[, 
            name := str_c(name, "_", seq_len(.N)),
            by=name
            ][]
}

#################################               ANALYSIS FUNCTIONS

binTableFromList <- function (input) {
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(magrittr)))
  elements <- unique(unlist(input))
  data <- lapply(input, function(x, e=elements){e %in% x %>% as.numeric}) %>% as.data.table()
  data[, Gene:=elements][]
}


phylopic_uid_item_safer <- function (name) {
  x <- gsub("[^a-zA-Z]+", "+", name)
  url <- paste0("http://phylopic.org/api/a/name/search?text=", 
                x, "&options=scientific+json")
  results <- jsonlite::fromJSON(url)$result
  if (length(results)==0){
    return(NA)
  } else{
    res <- results[[1]]
    for (id in res$uid) {
      uid <- ggimage:::phylopic_valid_id(id)
      if (!is.na(uid)) 
        break
    }
    return(uid)}
}

phylopic_uid_vector <- function (name) {
  res <- sapply(name, phylopic_uid_item_safer)
  return(res)
}


lineage.pathway.fdr.table <- function(path="output/ORA"){
  suppressMessages(suppressWarnings(require(tidyr)))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(readr)))
  suppressMessages(suppressWarnings(require(stringr)))
  options(readr.num_columns = 0)
  ORAs <- list.files(path, pattern="enrichment_results.*", recursive = T, full.names = T)
  tibble(
      Model = ORAs %>% sapply(str_extract, "(?<=ORA/).+?(?=/)"),
      RBHB = ORAs %>% sapply(str_extract, "RBB|evidenced"),
      is.dynamic = if_else(str_detect(Model, "_dyn"), "Yes", "No"),
      Database = ORAs %>% sapply(str_extract, "(?<=pathway_)[A-Za-z]+(?=-)"),
      FDR.Threshold = ORAs %>% sapply(str_extract, "(?<=FDR_)[0-9]\\.[0-9]+(?=/)"),
      Direction.lineage = ORAs %>% sapply(str_extract, "(?<=Project_).+?(?=/)"),
      Pathway.EnrichRatio.FDR.Genes.Size = ORAs %>% sapply(
        ., . %>% 
          read_tsv %>% 
          select(description, enrichmentRatio, FDR, userId, size) %>% 
          mutate_all(. %>% str_replace_all(",|;", ":")) %>% 
          unite("PEF", description, enrichmentRatio, FDR, userId, size, sep = ",") %>% 
          pull(PEF) %>% 
          paste0(collapse=";")),
      fname = ORAs
    ) %>% 
    separate(Direction.lineage, c("Direction", "Ancestor", "to", "Node"), sep="_", extra = "merge") %>% 
    separate_rows(Pathway.EnrichRatio.FDR.Genes.Size, sep=";") %>% 
    separate(Pathway.EnrichRatio.FDR.Genes.Size, c("Pathway","EnrichRatio","FDR", "Genes", "Size"), sep=",") %>% 
    group_by(Model, Database, Direction, Ancestor, Node, Pathway, EnrichRatio, Genes, RBHB) %>% 
    summarize(is.Dynamic = unique(is.dynamic) %>% paste0(collapse = ","), FDR.Threshold=min(FDR.Threshold), FDR=unique(FDR), Size=unique(Size) %>% paste(collapse=";")) %>% 
    ungroup %>% 
    mutate(Genes=str_replace_all(Genes, ":", ";")) %>% 
    arrange(Ancestor, Node, FDR, Pathway)
}

dynToCN <- function(val){
  val <- as.character(val)
  scale = c(NA, 1:36) %>% 
    set_names(c("-", as.character(0:9), LETTERS))
  #print(val)
  return(scale[[val]])
}

has.increased.dyn <- function(a, anc) {
  if (is.na(a) | is.na(anc)){
    return(FALSE)
  } else{
    return(a > anc)
  }
}

has.stable.dyn <- function(a, anc) {
  if (is.na(a) | is.na(anc)){
    return(FALSE)
  } else{
    return(a == anc)
  }
}

has.decreased.dyn <- function(a, anc) {
  if (is.na(a) | is.na(anc)){
    return(FALSE)
  } else{
    return(a < anc)
  }
}

genes.increased.dyn <- function(terminal, ancestor, gl){
  # print(nrow(terminal) == 1)
  # print(nrow(ancestor) == 1)
  gl %>% 
    lapply(
      ., 
      function(gene){
        # print(gene %in% names(terminal))
        # print(gene %in% names(ancestor))
        has.increased.dyn(a=terminal[[gene]], anc=ancestor[[gene]])
        }
    ) %>% 
    .[which(.==T)] %>% 
    names
}

getTPM <- function(genome){
  dir("../output/StringTie/final", pattern=genome, full.names = T) %>% 
    lapply(., . %>% read_tsv(col_names=F, comment="#") %>%
             # colnames)
             filter(str_detect(X9, "TPM")) %>%
             mutate(
               transcriptID = str_extract(X9, '(?<=transcript_id ")[A-Za-z0-9.]+(?=")'),
               TPM = str_extract(X9, '(?<=TPM ").+(?=")'),
               geneID = str_extract(X9, '(?<=gene_id ")[A-Za-z0-9.]+(?=")')
             )
    )%>%
    rbindlist %>%
    group_by(geneID) %>%
    summarize(TPM = max(as.numeric(TPM)))
}

getTPM.dt <- function(genome){
  suppressMessages(suppressWarnings(require(data.table)))
  dt <- dir("../output/StringTie/final", pattern=genome, full.names = T) %>% 
    lapply(
      ., 
      function(x){
        fread(x, header = F, select=c(3,9), col.names = c("feature", "attributes"))[
          feature=="transcript"
          ][,
            .(
              #transcriptID = str_extract(attributes, '(?<=transcript_id ")[A-Za-z0-9.]+(?=")'),
              TPM = str_extract(attributes, '(?<=TPM ").+(?=")'),
              geneID = str_extract(attributes, '(?<=gene_id ")[A-Za-z0-9.]+(?=")')
              )
          ]
        }
    ) %>%
    rbindlist()
  dt[, .(TPM = max(as.numeric(TPM))), by = geneID]
}
class.pathway <- function(Pathway){
  Pathway %>% 
    sapply(., function(P){
      ifelse(
        str_detect(P, "APC"),
        "APC-related Pathways",
        ifelse(
          str_detect(P, "[Pp]53"),
          "TP53 Pathway",
          ifelse(
            str_detect(P, "[Cc]heckpoint"),
            "Cell Cycle Checkpoint",
            ifelse(
              str_detect(P, "\\bG1(/S)?\\b|S Phase"),
              "Cell Cycle: G1, G1/S",
              ifelse(
                str_detect(P, "G2"),
                "Cell Cycle: G2",
                ifelse(
                  str_detect(P, "[Mm]itotic|M Phase|[Ss]ister [cC]hromatid"),
                  "Cell Cycle: Mitotic",
                  ifelse(
                    str_detect(P, "Cell [Cc]ycle|CDK"),
                    "Cell Cycle: General",
                    ifelse(
                      str_detect(P, "DNA.*[Dd]amage|[Rr]epair|\\b[Ii][Rr]\\b|Double Strand Break|Translesion synthesis|Homologous Recombination|Single Strand Annealing|Non.?[Hh]omologous [Ee]nd.[Jj]oining|NHEJ|NER\\b"),
                      "DNA Damage Response & Repair",
                      ifelse(
                        str_detect(P, "[Cc]ell [Dd]eath|Apoptosis|Ferroptosis|Necrosis"),
                        "Cell Death Pathways",
                        ifelse(
                          str_detect(P, "[Ss]enescence"),
                          "Senescence",
                          ifelse(
                            str_detect(P, "ATM|ATR|A[Kk][Tt]"),
                            "ATM, ATR, and AKT Pathways",
                            ifelse(
                              str_detect(P, "Pathways in Cancer|\\b[Rr][Aa][Ss]\\b|\\bRb\\b|Retinoblastoma|\\bPTEN\\b|Tumor suppressor|PI3"),
                              "Other Cancer-Related Pathways",
                              ifelse(
                                str_detect(P, "(?!Fluid)(\\b[Ss]tress\\b|\\b[Hh]ypoxia\\b|Oxidative damage|\\bUPR\\b)"),
                                "Cell Stress Response",
                                ifelse(
                                  str_detect(P, "TOR"),
                                  "mTOR Pathways",
                                  ifelse(
                                    str_detect(P, "Longevity|sirtuin|aging|Telomere"),
                                    "Longevity Regulating Pathway",
                                    ifelse(
                                      str_detect(P, "[Oo]ncogen"),
                                      "Oncogenic Pathway",
                                      ifelse(
                                        str_detect(P, "Insulin|IGF"),
                                        "Insulin- and IGF-Related Pathways",
                                        "Other"
                                      )
                                      )
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            ) 
          )
        )
      )})
}

graph.ORA.dotplot <- function(ORA.df, plotly=F){
  p <- ORA.df %>% 
    ggplot(aes(
      y=EnrichRatio, 
      x=FDR, 
      size=Size, 
      label = Pathway, 
      # text=GeneList, 
      color=Class
    )) + 
    geom_point() + #position = position_jitter(width = .02, height=.02, seed = NULL)) +
    # guides(col=guide_legend())+
    # facet_wrap(~Database) + 
    gghighlight(Class != "Other", keep_scales = T, calculate_per_facet = T, use_direct_label = F, label_key = NULL)+
    # scale_alpha_manual(values = c("Other"=0.20, 1), guide=F) +
    scale_x_log10(guide=guide_axis(title="FDR (log)")) +
    scale_y_log10(guide=guide_axis(title="Enrichment Ratio (log)")) +
    #geom_text_repel(data=. %>% filter(FDR<=1e-4|Class %in% c("TP53 Pathway", "Senescence", "Cell Death Pathways")), nudge_y = 0.5, show.legend = F, force = 15, size=3) +
    geom_text_repel(data=. %>% filter(FDR<=1e-4), nudge_y = 0.5,
                    show.legend = F, force = 15, size=3)+#, position = position_jitter(width = .02, height=.02, seed = NULL)) +
    geom_text_repel(data=. %>% filter(EnrichRatio>=7), nudge_y = -0.2, 
                    force = 15, show.legend=F, size=3)+#, position = position_jitter(width = .02, height=.02, seed = NULL)) +
    # geom_text_repel(data=. %>% filter(Class %in% c("TP53 Pathway", "Senescence", "Cell Death Pathways")), nudge_y = 0.5, force = 10, show.legend=F) + 
    scale_color_d3(palette="category20") +
    theme_pubclean() + 
    labs_pubr() + 
    theme(
      text = element_text(face = "plain", 
                          colour = "black", lineheight = 0.9, 
                          hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), 
                          debug = FALSE), 
      axis.text = element_text(size = 6, 
                               colour = "black", face = "bold"),
      axis.title = element_text(size = 8, 
                                colour = "black", face = "bold"),
      plot.title = element_text(size = 8, 
                                colour = "black", lineheight = 1, face = "bold"), 
      legend.text = element_text(size = 6, 
                                 face = "plain", colour = "black"),
      legend.position = "bottom"
    ) +
    guides(color=guide_legend(title="Cell Cycle Pathway", title.position = "bottom", title.hjust = 0.5), size = guide_legend(title="Size (log)", title.position = "bottom", title.hjust = 0.5)) 
  if(plotly) {p %<>% ggplotly}
  return(p)
}


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


not.me <- function(attr1, attr2, threshold=4){
  sapply(
    seq_along(attr1),
    function(n, a1=attr1, a2=attr2, t=threshold){
      ifelse(a1[n] >= t, a2[n], NA)
    }
  )
}
yes.me <- function(attr1, attr2, threshold=4){
  sapply(
    seq_along(attr1),
    function(n, a1=attr1, a2=attr2, t=threshold){
      ifelse(a1[n] < t, a2[n], NA)
    }
  )
}

colorado <- function(src, boulder) {                        # From https://stackoverflow.com/a/39695859, pun is too good to change
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- boulder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}

highlightTextSize <- function(src, boulder, sizes=c(4,8)) {
  if (!is.factor(src)) src <- factor(src)                 
  src_levels <- levels(src)                               
  brave <- boulder %in% src_levels                        
  if (all(brave)) {                                       
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels))
    b_vec <- rep(sizes[1], length(src_levels))            
    b_vec[b_pos] <- sizes[2]                              
    b_vec                                                   
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}

highlightTextColor <- function(src, boulder, color=c("grey","black")) {
  if (!is.factor(src)) src <- factor(src)                 
  src_levels <- levels(src)                               
  brave <- boulder %in% src_levels                        
  if (all(brave)) {                                       
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels))
    b_vec <- rep(color[1], length(src_levels))            
    b_vec[b_pos] <- color[2]                              
    b_vec                                                   
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}
# runGLS <- function(df, pcor = bm){
#   # print(class(df))
#   # print(head(df))
#   bla <- gls(
#     media.CN ~ Genome + N50 + meanSequenceLength + Completeness + Total+ pc.Ortho + Average + , 
#     data = df, 
#     correlation=pcor
#   ) 
#   # print(bla)
#   # print(class(bla))
#   return(bla)
# }
tidyGLS <- function(df){
  # print(names(df))
  # print(class(df$results))
  a <- tidy(df, conf.int=T)
  print(a)
  return(a)
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

text.ratio.nodeByMRCA.ancNode <- function(nodeA, outNode, tree_table=tree.Atlantogenata.table, nodeid=T){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  if (nodeid){
    nodeA = nodeid(tree_table, nodeA)
    outNode = nodeid(tree_table, outNode)
  }
  real.node = tree_table %>% filter(node == MRCA(tree_table, nodeA, outNode) %>% pull(node)) %>% pull(node)
  anc.node = tree_table %>% filter(node == MRCA(tree_table, nodeA, outNode) %>% pull(node)) %>% pull(parent)
  bsize.A <- tree_table %>% 
    filter(node == real.node) %>%
    pull(Median)
  # print(bsize.A)
  bsize.B <- tree_table %>% 
    filter(node == anc.node) %>% 
    pull(Median)
  # print(bsize.B)
  ratio.bsize <- bsize.A/bsize.B
  if (!is.numeric(ratio.bsize)) {warning(sprintf("%s is not numeric", ratio.bsize))}
  ratio.bsize <- ifelse(
    ratio.bsize<1, 
    divide_by(1, ratio.bsize) %>% round(digits=2) %>% paste0(., "x decrease"), 
    ratio.bsize %>% round(digits=2) %>% paste0(., "x increase")
  )
  return(ratio.bsize)
}

text.ratio.node.ancnode.byMRCA <- function(nodeA, outNode, tree_table=tree.Atlantogenata.table, nodeid=T){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  if (nodeid){
    nodeA = nodeid(tree_table, nodeA)
    outNode = nodeid(tree_table, outNode)
  }
  anc.node = MRCA(tree_table, nodeA, outNode) %>% pull(node)
  bsize.A <- tree_table %>% 
    filter(node == nodeA) %>%
    pull(Median)
  # print(bsize.A)
  bsize.B <- tree_table %>% 
    filter(node == anc.node) %>% 
    pull(Median)
  # print(bsize.B)
  ratio.bsize <- bsize.A/bsize.B
  if (!is.numeric(ratio.bsize)) {warning(sprintf("%s is not numeric", ratio.bsize))}
  ratio.bsize <- ifelse(
    ratio.bsize<1, 
    divide_by(1, ratio.bsize) %>% round(digits=2) %>% paste0(., "x decrease"), 
    ratio.bsize %>% round(digits=2) %>% paste0(., "x increase")
  )
  return(ratio.bsize)
}

text.ratio.node.ancnode.byNodes <- function(nodeA, anc.node, tree_table=tree.Atlantogenata.table){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  bsize.A <- tree_table %>% 
    filter(node == nodeA) %>%  
    pull(Median)
  # print(bsize.A)
  bsize.B <- tree_table %>% 
    filter(node == anc.node) %>% 
    pull(Median)
  # print(bsize.B)
  ratio.bsize <- bsize.A/bsize.B
  if (!is.numeric(ratio.bsize)) {warning(sprintf("%s is not numeric", ratio.bsize))}
  ratio.bsize <- ifelse(
    ratio.bsize<1, 
    divide_by(1, ratio.bsize) %>% round(digits=2) %>% paste0(., "x decrease"), 
    ratio.bsize %>% round(digits=2) %>% paste0(., "x increase")
  )
  return(ratio.bsize)
}

text.ratio.node.ancnode.byLabel <- function(node.label, ancNode.label, tree_table=tree.Atlantogenata.table){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  bsize.A <- tree_table %>%
    mutate(label=str_replace(label, " ", ".")) %>%  
    filter(label==node.label) %>%
    pull(Median)
  # print(bsize.A)
  bsize.B <- tree_table %>% 
    mutate(label=str_replace(label, " ", ".")) %>% 
    filter(label==ancNode.label) %>% 
    pull(Median)
  # print(bsize.B)
  ratio.bsize <- bsize.A/bsize.B
  if (!is.numeric(ratio.bsize)) {warning(sprintf("%s is not numeric", ratio.bsize))}
  ratio.bsize <- ifelse(
    ratio.bsize<1, 
    1/ratio.bsize %>% round(digits=2) %>% paste0(., "x decrease"), 
    ratio.bsize %>% round(digits=2) %>% paste0(., "x increase")
  )
  return(ratio.bsize)
}

text.ratio.node.ancnode.atLabel <- function(node.label, tree_table=tree.Atlantogenata.table){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  nodeA <- nodeid(tree_table, node.label)
  anc.node <- parent(tree_table, node.label) %>% pull(node)
  bsize.A <- tree_table %>%
    filter(node==nodeA) %>%
    pull(Median)
  # print(bsize.A)
  bsize.B <- tree_table %>% 
    filter(node==anc.node) %>% 
    pull(Median)
  # print(bsize.B)
  ratio.bsize <- bsize.A/bsize.B
  if (!is.numeric(ratio.bsize)) {warning(sprintf("%s is not numeric", ratio.bsize))}
  ratio.bsize <- ifelse(
    ratio.bsize<1, 
    1/ratio.bsize %>% round(digits=2) %>% paste0(., "x decrease"), 
    ratio.bsize %>% round(digits=2) %>% paste0(., "x increase")
  )
  return(ratio.bsize)
}


text.bodysize.ancNode.byMRCA <- function(nodeLabelA, nodeLabelB, tree_table=tree.Atlantogenata.table){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  anc.node = MRCA(tree_table, nodeLabelA, nodeLabelB) %>% pull(node)
  bsize <-  tree_table %>% 
    filter(node == anc.node) %>% 
    pull(Median) %>% 
    divide_by(1000)
  # print(bsize)
  # print(class(bsize))
  bsize <- ifelse(
    bsize < 1,
    bsize %>% multiply_by(1000) %>% round(., digits = 2) %>% paste0(., "g"),
    bsize %>% round(digits = 2) %>% paste0(., "kg")
  )
    return(bsize)
}

text.bodysize.ancNode.byMRCA.CI <- function(nodeA, nodeB, tree_table=tree.Atlantogenata.table){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  anc.node = MRCA(tree_table, nodeA, nodeB) %>% pull(node)
  bsize <-  tree_table %>% 
    filter(node == anc.node) %>% 
    select(CI_low=`95%CI_Low`, CI_High=`95%CI_High`) %>% 
    transmute(CI_range = paste0(CI_low %>% exp %>% divide_by(1000) %>% round(digits=2), "kg-", CI_High %>% exp %>% divide_by(1000) %>% round(2), "kg")) %>% 
    pull(CI_range)
  return(bsize)
}

text.bodysizeCI.node.byLabel <- function(node.label, tree_table=tree.Atlantogenata.table){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  bsize <-  tree_table %>% 
    filter(label == node.label) %>% 
    select(CI_low=`95%CI_Low`, CI_High=`95%CI_High`) %>% 
    transmute(CI_range = paste0(CI_low %>% exp %>% divide_by(1000) %>% round(digits=2), "kg-", CI_High %>% exp %>% divide_by(1000) %>% round(2), "kg")) %>% 
    pull(CI_range)
  return(bsize)
}

text.bodysize.byLabel <- function(node.label, tree_table=tree.Atlantogenata.table){
  tree_table %<>% mutate(label=str_replace(label, " ", "."), Median = Median %>% exp)
  bsize <-  tree_table %>% 
    filter(label==node.label) %>% 
    pull(Median) %>% 
    divide_by(1000) %>% 
    round(digits=2) %>% 
    paste0(., "kg")
  return(bsize)
}

as.percentage <- function(decimal){
  scales::label_percent(accuracy = 0.01)(decimal)
}
subtract.rows <- function(r1, r2, exclude_vars){
  keep.label <- r1 %>% select(parent, node, branch.length, label, starts_with(exclude_vars))
  subtractand <- r1 %>% select(-parent, -node, -branch.length, -label, -starts_with(exclude_vars))
  subtractor <- r2 %>% select(-parent, -node, -branch.length, -label, -starts_with(exclude_vars))
  new.vals <- subtractand-subtractor
  new.row <- cbind(keep.label, new.vals)
  return(new.row)
}

gene.increase.tree <- function(tree.data, node=NULL, root.node=NULL, parent=NULL, exclude="lnSize"){
  # Lol look at me programming a tree traversal algorithm
  # first get root node
  if (is.null(root.node)){
    final.table <- tree.data %>% select(parent, node, branch.length, label)
    root.node = rootnode(tree.data) %>% pull(node) %>% unique
    my.node = root.node
  } else{
    my.node = node
  }
  parent.node = tree.data %>% filter(node==my.node) %>% pull(parent) %>% unique
  child.nodes = child(tree.data, my.node) %>% pull(node)
  if(length(child.nodes) != 0){
    child.rows <- lapply(
      child.nodes,
      function(n, t = tree.data, r=root.node){
        gene.increase.tree(tree.data = t, node=n, root.node = r)
      }
    ) %>% bind_rows() 
  }else{
    child.rows = tree.data %>% filter(label=="There's no way you'll match this as a label unless you're trying really hard.")
  }
  if(my.node != root.node){
    t1 = tree.data %>% filter(node==my.node)
    t2 = tree.data %>% filter(node==parent.node)
    my.row <- subtract.rows(t1, t2, exclude_vars = exclude)
  }else{
    t1 = tree.data %>% filter(node==my.node)
    t2 = tree.data %>% filter(node==root.node)
    my.row <- subtract.rows(t1, t2, exclude_vars = exclude)
  }
  new.table <- bind_rows(my.row, child.rows) %>% as_tibble
  if (my.node == root.node){
    new.table <- new.table %>% as_tibble
    new.table <- full_join(final.table, new.table, by=c("parent", "node", "branch.length", "label"))
  }
  return(new.table)
}

branch.size.ratio <- function(n, tree.data){
  if(!is_tibble(tree.data)){tree.data %<>% as_tibble}
  node.size <- tree.data %>% filter(node==n) %>% pull(lnSize) %>% unique() %>% exp
  anc.node <- tree.data %>% filter(node==n) %>% pull(parent) %>% unique()
  # print(label)
  # print(anc.node)
  anc.node.size <- tree.data %>% filter(node==anc.node) %>% pull(lnSize) %>% unique() %>% exp
  ratio <- node.size %>% divide_by(anc.node.size) %>% log %>% round(digits=2)
  return(ratio)
}

get.ancestor.label <- function(n, tree.data){
  if(!is_tibble(tree.data)){tree.data %<>% as_tibble}
  anc.node <- tree.data %>% filter(node==n) %>% pull(parent) %>% unique()
  ancestor.label <- tree.data %>% filter(node==anc.node) %>% pull(label)
  return(ancestor.label)
}
