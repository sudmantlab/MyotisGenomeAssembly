Rules
=====


This folder contains the Snakemake rules used during assembly and QC. 

Generally, each Snakefile contains rules pertaining to a certain aspect of the assembly. 

- `Hifi.smk` contains rules associated with generating and assembling HiFi reads;
- `genomeQC.smk` contains rules associated with generating QC metrics at each step.
