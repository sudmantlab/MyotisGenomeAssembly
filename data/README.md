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

- `lifehistory/`: Contains data sourced from different databases used for longevity and other
life history analyses
  - `AnAge/`: Dataset containing an assortment of life history data, centered on maximum lifespan.
  - `PanTHERIA/`: A species-level database of life history, ecology, and geography of extant and
recently extinct mammals.

- `GenomeQCStats/`: Contains the final genome stats and quality control metrics, generated in `analyses/QC`.
  - `assembly_stats/`: Contains cleaned output table of genome summary stats from Peter's
`assembly_stats` script
    - `2022-11-4/`
      - `chiroptera.assembly_stats.tsv`: Table of genome summary statistics for our genomes, plus
other published bat and select other genomes.

- `USGS_SpeciesRanges/`: Contains species ranges pulled from the
[USGS GAP Analysis Project](https://www.usgs.gov/programs/gap-analysis-project), sorted by species.

- `tree/`: Contains trees used in this project, either published or generated in `analyses/phylogeny`.
  - `species_timetree.nwk`: Original species tree pulled from [timetree.org](https://timetree.org)
  - 'Myotis_orthofinder.nwk`: Tree generated from orthofinder
