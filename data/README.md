Data folder for analyses
========================

This folder contains all the FINAL data sources that we will be using for different analyses. 

Folders should be organized by the analysis/program/purpose of the data type, followed by date, 
ending in the final files. When adding folders here, please update the "CURRENT 
STRUCTURE" section below to include a description of the folders' contents.

*Important* - if you're adding a published dataset here, also make sure to include a `.bib` file 
with the citation for the data set!

CURRENT STRUCTURE
-----------------

- `lifehistory/`: contains data sourced from different databases used for longevity and other 
life history analyses
  - `AnAge/`: Dataset containing an assortment of life history data, centered on maximum lifespan.
  - `PanTHERIA/`: A species-level database of life history, ecology, and geography of extant and 
recently extinct mammals.

- `GenomeQCStats/`: contains the final genome stats and quality control metrics, generated in `analyses/QC`.
  - `assembly_stats/`: contains cleaned output table of genome summary stats from Peter's 
`assembly_stats` script
    - `2022-11-4/`
      - `chiroptera.assembly_stats.tsv`: Table of genome summary statistics for our genomes, plus 
other published bat and select other genomes.



