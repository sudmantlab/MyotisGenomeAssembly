# Bat Genome Analyses

## Folder Structure

- data

  All the raw and cleaned data should be placed and organized in the `data` folder. 
  Currently this includes the genomes located in `data/genomes`, the transcriptomes 
  in `data/trinity`, and mapped reads in the `data/[data_type]-bam` folders.

- analysis

  Each self-contained analysis folder should be uploaded into this folder. Ideally,
  all code and paths for inputs, etc. should be relative within this folder, and  
  a data folder inside the analysis folder should symlink to data in the main data 
  folder for consistency and to save space. Any finalized outputs for different 
  analyses should be symlinked to the main data folder as well.

- code_bib

  As we use code that's either published or on github that's citable, it might be 
  helpful to have a place to store all of the `.bib` files for later citing. 
  
- envs

  I use/abuse conda for all of my analyses, as I assume others do, and so having a 
  folder to share conda environments will be really useful for reproducibility.



