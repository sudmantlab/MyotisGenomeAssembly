if (!requireNamespace("WebGestaltR", quietly = TRUE))
  install.packages("WebGestaltR", repos='http://cran.us.r-project.org')
suppressMessages(suppressWarnings(require(WebGestaltR, quietly=T)))


run_ORA <- function(geneSet, refGeneSet, outputDirectory, enrichmentDB, fdrThr, projName){
  WebGestaltR(
    enrichMethod = "ORA",
    enrichDatabase = enrichmentDB,
    organism = "hsapiens",
    interestGeneFile = geneSet,
    interestGeneType = "genesymbol",
    referenceGeneFile = refGeneSet,
    referenceGeneType = "genesymbol",
    maxNum=1000,
    sigMethod="top",
    topThr=100,
    fdrThr=as.numeric(fdrThr),
    reportNum=100,
    isOutput = T,
    outputDirectory = outputDirectory,
    hostName="http://www.webgestalt.org/",
    projectName=projName
  )
}

dir.create(snakemake@output[["outputDirectory"]], recursive = T, showWarnings = T)

run_ORA(geneSet=snakemake@input[["geneSet"]], refGeneSet=snakemake@input[["refGeneSet"]], outputDirectory=snakemake@params[["rootDir"]], 
        enrichmentDB=snakemake@wildcards[["enrichmentDB"]], fdrThr = snakemake@wildcards[["fdr"]], projName = snakemake@wildcards[["list"]])
