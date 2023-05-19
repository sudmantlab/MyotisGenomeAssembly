library(rmarkdown)
library(rticles)
render(snakemake@input[[1]], output_format="all", envir=new.env())
