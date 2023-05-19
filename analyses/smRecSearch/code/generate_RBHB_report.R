#!/usr/bin/env Rscript

library(knitr)
library(rmarkdown)
library(optparse)

option_list = list(
  make_option(c("-s", "--species"), type="character", default="Mus musculus",
              help="Species Name  [eg: %default]", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="Lab Mouse",
              help="Common Name [eg: %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="mm10",
              help="Genome [eg: %default]", metavar="character"),
  make_option("--template", type="character", default="./analysis/template_RBHB_analysis.Rmd",
              help="Template file used for report [default: %default]", metavar="character"),
  make_option("--outDir", type="character", default="./docs",
              help="Output Directory [default: %default]", metavar="character"),
  make_option("--workDir", type="character", default="./",
              help="Working Directory [default: %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

opt$template = normalizePath(opt$template)
opt$outDir = normalizePath(opt$outDir)

setwd(opt$workDir)

render(
  opt$template,
  params = list(
    name = opt$name,
    species = opt$species,
    genome = opt$genome
    ),
  output_file = paste0(opt$outDir, "/Duplicate_Analysis_", opt$genome, ".html")
  )

