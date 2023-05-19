suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(magrittr)))
suppressMessages(suppressWarnings(require(stringr)))
suppressMessages(suppressWarnings(require(readr)))

files <- snakemake@input[["geneCopyTables"]]
genomeTable <- fread(snakemake@input[["genomeTable"]])
genomeTable[, Genome:=BestGenome]
setkey(genomeTable, Genome)

genomes <- sapply(files, basename) %>% sapply(str_remove, "_geneCopyTable_(RBB|evidenced).*.tsv$")

names(files) <- genomes

getGeneTable <- function(file_path, genome){
  dt <- fread(file_path)#[, .(Gene, CopyECNC=do.call(paste,c(.SD, sep=","))), .SDcols=c("Copy", "ECNC")]
  # print(nrow(dt))
  # print(names(dt))
  # print(data.table::transpose)
  # print(packageVersion("data.table"))
  #setnames(dt, paste0(genome, "_", names(dt)))
  #dt <- data.table::transpose(l = dt, make.names = "Gene")
  dt[,Genome:=genome]
  setkey(dt, Genome)
  setcolorder(dt)[]
}

DT.replace.na.qmark.DT = function(DT) {
  for (j in seq_len(ncol(DT))){set(DT,which(is.na(DT[[j]])),j,"?")}
}

geneCopyTable <- lapply(
  genomes, 
  function(genome, filepath_list = files){
    # print(genome)
    getGeneTable(filepath_list[[genome]], genome)
  }
) %>%
  rbindlist

setkey(geneCopyTable, Genome)

code.duplicate <- function(Copy, ECNC){
  ifelse(
    Copy==0,
    return(0),
    ifelse(
      Copy == 1,
      return(1),
      ifelse(
        ECNC >= 1.5,
        return(2),
        return(1)
      )
    )
  )
}

code.duplicate.dynamic <- function(Copy, ECNC){
  # Round ECNC up
  ECNC = ECNC %>% sum(., ifelse(.>round(.), 0.01, 0)) %>% round()
  # set Copy Number to the smallest of the two
  CN = ifelse(ECNC<Copy, ECNC, Copy) - 1
  # CORRECTION: IQTree does NOT allow >32 states!!!
  result <- ifelse(
    is.na(CN),
    "?",
      ifelse(
        CN < 10, 
        as.character(CN), 
        ifelse(
          CN < 31,
          LETTERS[CN-9],
          LETTERS[31-9]
        )
      )
    )
  if(is.na(result)){
    print(CN)
  } else {return(result)}
}

fwrite(x = geneCopyTable, file=snakemake@output[["longTable"]])


geneCopyTable.wide <- merge(
  dcast(geneCopyTable[,.(value = code.duplicate(Copy,ECNC) %>% as.character()), by=.(Genome,Gene)], Genome ~ Gene),
  genomeTable[, strsplit(as.character(Genome), ",", fixed=TRUE), by = .(label, Genome)][,.(Genome = V1, label)]
)


geneCopyTable.wide[, Genome:=NULL]
setkey(geneCopyTable.wide, label)
setcolorder(geneCopyTable.wide)

fwrite(x = geneCopyTable.wide, file=snakemake@output[["wideTable"]])


# Should always return True:
# length(geneCopyTable.wide[, .(label=label, k = do.call(paste, c(.SD, sep = ""))), .SDcols = -"label"][,str_length(k)]) == 1
DT.replace.na.qmark.DT(geneCopyTable.wide)

first_line = sprintf("%s %s", nrow(geneCopyTable.wide), ncol(geneCopyTable.wide)-1)

c(
  first_line,
  geneCopyTable.wide[, .(label=label, k = do.call(paste, c(.SD, sep = ""))), .SDcols = -"label"][, do.call(paste, c(.SD, sep = " "))]
) %>% 
  write_lines(snakemake@output[["codedphylip"]])



geneCopyTable.wide.dyn <- merge(
  dcast(geneCopyTable[,.(value = code.duplicate.dynamic(Copy,ECNC) %>% as.character()), by=.(Genome,Gene)], Genome ~ Gene),
  genomeTable[, strsplit(as.character(Genome), ",", fixed=TRUE), by = .(label, Genome)][,.(Genome = V1, label)]
)

geneCopyTable.wide.dyn[, Genome:=NULL]
setkey(geneCopyTable.wide.dyn, label)
setcolorder(geneCopyTable.wide.dyn)

fwrite(x = geneCopyTable.wide.dyn, file=snakemake@output[["wideTable_dyn"]])

first_line.dyn = sprintf("%s %s", nrow(geneCopyTable.wide.dyn), ncol(geneCopyTable.wide.dyn)-1)

# Should always return True:
# length(unique(geneCopyTable.wide.dyn[, .(label=label, k = do.call(paste, c(.SD, sep = ""))), .SDcols = -"label"][,str_length(k)])) == 1

c(
  first_line.dyn,
  geneCopyTable.wide.dyn[, .(label=label, k = do.call(paste, c(.SD, sep = ""))), .SDcols = -"label"][, do.call(paste, c(.SD, sep = " "))]
) %>% 
  write_lines(snakemake@output[["codedphylip_dyn"]])


colnames(geneCopyTable.wide[, !"label"]) %>% 
  write_lines(snakemake@output[["phylipGeneList"]])

save.image("debug_consolidateGeneCopyTable.RData")