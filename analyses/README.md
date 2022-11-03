Bat Genome Analyses
===================

Folder Structure
------------------

- `data`

  All the raw and final/clean data should be symlinked or copied into the `data` folder. 
  Currently this includes the genomes located in `data/genomes`, the transcriptomes 
  in `data/trinity`, and mapped reads in the `data/[data_type]-bam` folders.

- `analysis`

  Each self-contained analysis folder should be uploaded into this folder. Ideally,
  all code and paths for inputs, etc. should be relative within this folder, and  
  a data folder inside the analysis folder should symlink to data in the main data 
  folder for consistency and to save space. Any finalized outputs for different 
  analyses should be symlinked to the main data folder as well.
  **Sometimes rather than a symlink, an analysis folder will have its own 

- `code_bib`

  As we use code that's either published or on github that's citable, it might be 
  helpful to have a place to store all of the `.bib` files for later citing. 
  
- `envs` and `rules`

  I use/abuse conda for all of my analyses, as I assume others do, and so having a 
  folder to share conda environments will be really useful for reproducibility.
  
- `figures`
  
  A central place to hold all putatively finalized figures. 
  Internal folder structure should reflect the analyses folders, further organized by
  date of creation of the figures.

Analyses
--------

- annotation  ![nearestNeighborTE](https://img.shields.io/badge/Status-Ongoing-yellow)
  
  Genome annotation using Funannotate and TOGA.

- cactus  ![nearestNeighborTE](https://img.shields.io/badge/Status-Pending-blue)
  
  Cactus-based pangenome alignment using `progressiveCactus` and `minigraph-cactus`.

- COSMIC  ![nearestNeighborTE](https://img.shields.io/badge/Status-Pending-blue)
  
  Analyses of selection and coding changes specifically at COSMIC genes.

- heterozygosity  ![nearestNeighborTE](https://img.shields.io/badge/Status-Verify-red)
  
  Calling heterozygosity and identifying ROHs across each genome using a variety of methods.

- minimap2  ![nearestNeighborTE](https://img.shields.io/badge/Status-Done-green)
  
  Generates BAM files of PacBio Reads to each genome.
  
- nearestNeighborTE	![nearestNeighborTE](https://img.shields.io/badge/Status-Pending-blue)

  A simple pipeline to look for overrepresented co-localization of TEs and genes.

- pangenome  ![nearestNeighborTE](https://img.shields.io/badge/Status-Pending-blue)

  Contains all code, etc necessary to build the pan-Myotis genome graphs and outputs.

- QC  ![nearestNeighborTE](https://img.shields.io/badge/Status-Ongoing-yellow)
  
  Contains all the assembly statistics for genomes used in this paper.

- smRecSearch  ![nearestNeighborTE](https://img.shields.io/badge/Status-Verify-red)
  
  Running my Reciprocal Best-Hit BLAT pipeline to have a quick-and-dirty estimation of gene copy numbers.
