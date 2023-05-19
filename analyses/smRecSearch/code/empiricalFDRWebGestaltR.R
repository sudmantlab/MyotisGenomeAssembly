library(tidyverse)
options(readr.num_columns = 0)
library(data.table)
library(magrittr)
library(WebGestaltR)
library(parallel)
library(ggpubr)
source("../code/generalFunctions.R")
source("../code/R_rainclouds.R")
source("../code/summarySE.R")

dir.create("../output/empiricalORA", recursive = T)

no_cores <- detectCores() - 2 # I want to still use my laptop...
cl <- makeCluster(no_cores, outfile="")
clusterExport(cl, "WebGestaltR")

set.seed(1234)

reactome.ORA <- function(geneList, refList, is_output=F, fdr.thres=1){
  suppressWarnings(require(WebGestaltR, quietly = T))
  suppressMessages(suppressWarnings(WebGestaltR(
    enrichMethod = "ORA",
    enrichDatabase = "pathway_Reactome",
    organism = "hsapiens",
    interestGene = geneList,
    interestGeneType = "genesymbol",
    referenceGene = refList,
    referenceGeneType = "genesymbol",
    maxNum=1000,
    fdrMethod="BH",
    sigMethod="fdr",
    fdrThr = fdr.thres,
    reportNum=100,
    isOutput = is_output,
    # outputDirectory = outputDirectory,
    hostName="http://www.webgestalt.org/"
  )))
}

gene_list <- read_lines("../output/recBlastDBPrep/AvA_geneList_filtered.txt")
geneCopyTable <- read_csv("../output/geneCopyTable/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata-GeneCopyTable_RBB_filtered_long.csv")
genomeTable <- read_csv("../output/other/genomeTable.csv")

geneCopyTable.stats <- geneCopyTable %>% 
  group_by(Genome) %>% 
  summarize(N.Hits = sum(!is.na(ECNC) & !is.na(Copy)),
            N.Dup = sum(ECNC>=1.5 & Copy>1, na.rm = T))
set.sizes <-c(1,3,5,9,12,15,28,35,52,61,76,83,94,157,383,484,504,587,1000,1591, 2103, 3000, 4000, 5000) %>% 
  set_names(., as.character(.))

fdr.list <- c(
  # 0.1, 
  0.25, 
  0.5, 
  0.75, 
  1
  )

getGeneSet <- function(set.size, ref_list){
  sample(
    x = ref_list,
    size= set.size,
    replace=F
    )
}

clusterExport(cl, "getGeneSet")
clusterExport(cl, "reactome.ORA")
clusterExport(cl, "gene_list")
clusterExport(cl, "tibble")
clusterExport(cl, "%>%")
clusterExport(cl, "mutate")
clusterExport(cl, "select")
clusterExport(cl, "as.data.table")

run.ORAs <- function(sets_, run, ref_list_, fdrThreshold){
  parLapply(
  # lapply(
    cl,
    sets_, 
    function(x, y=ref_list_, r = run){
      gSet <- getGeneSet(x, y)
      ORA <- tryCatch(
        reactome.ORA(geneList = gSet, refList = y, fdr.thres = fdrThreshold),
        error = function(err){
          print(err)
          tibble(geneSet="", description=NA, link=NA, size=NA, overlap=NA, expect=NA, enrichmentRatio=NA, pValue=NA, FDR=NA, overlapId=NA, userId=NA)
          }
      )
      if (is.null(ORA)){
        print("ORA is null")
        ORA <- tibble(geneSet="", description=NA, link=NA, size=NA, overlap=NA, expect=NA, enrichmentRatio=NA, pValue=NA, FDR=NA, overlapId=NA, userId=NA)
      }
      ORA %>%
        mutate(N.set = x, Run = r) %>% 
        select(Run, N.set, description, pValue) %>% 
        as.data.table() %>% 
        return(.)
    }
  ) %>% 
    rbindlist()
}


replicate.run.ORAs <- function(fdr.threshold, n=1:5000, sets = set.sizes, ref_list = gene_list){
  lapply(
    n,
    function(x, a=sets, b=ref_list, c=fdr.threshold){
      run.ORAs(x, sets_ = a, ref_list_ = b, fdrThreshold = c)
      }
  ) %>% 
    rbindlist()
}

digest.data <- function(dt, fdr=1){
  empiricalFreq.dt <- dt[, 
                         .(Runs = max(Run), description=description, N.set=N.set)
                         ][, 
                           .(Overall.Freq = .N/max(Runs)), by = .(description, N.set)
                           ]
  setkey(empiricalFreq.dt, N.set, description)
  
  empiricalFreq2.dt <- dt[,
                          .(Run, N.set, description, pValue)
                          ][, 
                            .(N.Pathways = .N, description = description, pValue=pValue), by=.(Run, N.set)
                            ][,
                              .(freq.Description = .N/N.Pathways, pValue=pValue), by=.(Run, N.set, description)
                              ][,
                                .(pValue.mean = mean(pValue, na.rm = T), pValue.Median = median(pValue, na.rm = T), 
                                  pValue.sd = sd(pValue, na.rm = T),
                                  PerRun.freq.mean = mean(freq.Description, na.rm = T), PerRun.freq.Median = median(freq.Description, na.rm = T), 
                                  PerRun.freq.sd = sd(freq.Description, na.rm = T)), 
                                by = .(description, N.set)
                                ]
  setkey(empiricalFreq2.dt, N.set, description)
  
  empiricalFreq3.dt <- merge(empiricalFreq.dt, empiricalFreq2.dt)
  
  dt %>% fwrite(file = sprintf("../output/empiricalORA/empiricalORA_seed1234_FDR%s.csv", fdr*100))
  empiricalFreq3.dt %>% fwrite(file=sprintf("../output/empiricalORA/stats_empiricalORA_seed1234_FDR%s.csv", fdr*100))
  return(empiricalFreq3.dt)
}

for (fdr.level in fdr.list){
  print(fdr.level)
  clusterExport(cl, "fdr.level")
  data.set <- replicate.run.ORAs(fdr.threshold = fdr.level)
  empirical.ds <- digest.data(data.set, fdr = fdr.level)
  stats.empirical.ds <- empirical.ds[
    !is.na(pValue.mean), 
    .(N = .N, pval.mean=mean(pValue.mean, na.rm = T), pval.sd=sd(pValue.mean, na.rm = T), pval.median=median(pValue.mean, na.rm = T),
      perRun.mean=mean(PerRun.freq.mean, na.rm = T), perRun.sd=sd(PerRun.freq.mean, na.rm = T), perRun.median=median(PerRun.freq.mean, na.rm = T)), 
    by="N.set"
    ][, `:=`(pval.ci=qt(0.95,df=N-1)*pval.sd/sqrt(N), perRun.ci=qt(0.95,df=N-1)*perRun.sd/sqrt(N)), by="N.set"][]
  
  stats.empirical.ds %>% 
    fwrite(sprintf("../output/empiricalORA/pValueStats_empiricalORA_seed1234_FDR%s.csv", fdr.level))
  
  p <- empirical.ds %>%
  filter(N.set %in% c(1, 5, 25, 50, 100, seq.int(from=500, to=5000, by = 500))) %>%
  ggplot(
    aes(
      x=factor(N.set, levels = N.set %>% sort %>% unique), 
      y=PerRun.freq.mean
      )
    ) +
    geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2)+
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_point(
      data = stats.empirical.ds %>% filter(N.set %in% c(1, 5, 25, 50, 100, seq.int(from=500, to=5000, by = 500))), 
      aes(
        x = factor(N.set, levels = N.set %>% sort %>% unique), 
        y = perRun.mean
      ), 
      position = position_nudge(.25), 
      colour = "BLACK", 
      inherit.aes = F
      )+
    geom_errorbar(
      data = stats.empirical.ds %>% filter(N.set %in% c(1, 5, 25, 50, 100, seq.int(from=500, to=5000, by = 500))), 
      aes(
        x = factor(N.set, levels = N.set %>% sort %>% unique), 
        y = perRun.mean, 
        ymin = perRun.mean - perRun.ci, 
        ymax = perRun.mean + perRun.ci
        ), 
      position = position_nudge(.25), 
      colour = "BLACK", 
      width = 0.1, 
      size = 0.8, inherit.aes = F
      )+
    theme_pubclean()+
  scale_y_log10() +
  labs(x= "Number of Genes Sampled",
       y= "log(Frequency per Run of Pathway)",
       caption=sprintf("FDR: %s", fdr.level))
  p %>% ggsave(sprintf("../output/empiricalORA/PathwayFreq_empiricalORA_seed1234_FDR%s.svg", fdr.level), plot = ., device = "svg", width = 16, height = 9)

p <- empirical.ds %>%
  filter(N.set %in% c(1, 5, 25, 50, 100, seq.int(from=500, to=5000, by = 500))) %>%
  ggplot(aes(x=factor(N.set, levels = N.set %>% sort %>% unique), y=pValue.mean)) +
    geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2)+
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_point(
      data = stats.empirical.ds %>% filter(N.set %in% c(1, 5, 25, 50, 100, seq.int(from=500, to=5000, by = 500))), 
      aes(
        x = factor(N.set, levels = N.set %>% sort %>% unique), 
        y = pval.mean
        ), 
      position = position_nudge(.25), 
      colour = "BLACK", 
      inherit.aes = F
      )+
    geom_errorbar(
      data = stats.empirical.ds %>% filter(N.set %in% c(1, 5, 25, 50, 100, seq.int(from=500, to=5000, by = 500))), 
      aes(
        x = factor(N.set, levels = N.set %>% sort %>% unique), 
        y = pval.mean, 
        ymin = pval.mean-pval.ci, 
        ymax = pval.mean+pval.ci
        ), 
      position = position_nudge(.25), 
      colour = "BLACK", 
      width = 0.1, 
      size = 0.8, 
      inherit.aes = F
      )+
    theme_pubclean()+
  scale_y_log10() +
  labs(x= "Number of Genes Sampled",
       y= "log(P-Values)",
       caption=sprintf("FDR: %s", fdr.level))
  p %>% ggsave(sprintf("../output/empiricalORA/PValueDist_empiricalORA_seed1234_FDR%s.svg", fdr.level), plot = ., device = "svg", width = 16, height = 9)

}

stopCluster(cl)
