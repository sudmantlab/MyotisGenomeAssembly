Assembly of Chromosome-Level Genomes
====================================

Our genomes were assembled using a combination


Repo Structure: 
-------

- `data/`: Contains all initial data for assembly, including sequencing data. 
Also contains folders with symlinks to final products;
- `envs/`: Contains the Conda environment files defining the tools for assembly. 
Env files with dates are the explicit environments used in generating the assemblies;
- `output/`: Contains all the intermediate and final outputs of assembly;
- `rules/`: Contains the Snakefiles with rules for assembling and QC-ing the genomes;
- `Snakefile`: The master Snakefile for assembling the genomes;
- `pepsamples.tsv`: Table with all the sequencing metadata (derived from the files themselves);
- `
