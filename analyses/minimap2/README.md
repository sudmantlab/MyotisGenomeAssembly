This directory contains all of the PacBio-mapped BAM files to the final genome files. 
These BAMs will be used to call variants and assess gaps and/or assembly errors later on.

Directory Structure
===================

- `data/`
  - `HiFi-adapterFiltered [symlink]`: PacBio HiFi CCS BAM/FastQ files with SMRTBELL adapters removed (usually unnecessary, 100% a preemptive precaution)
  - `genomes [symlink]`: contains FASTA, index, and genome files for all genomes. 
  NCBI genomes have species names, while our genomes are named `mMyo[Spe]1`. 
  Note that currently as of 2022-09-08, `mMyo[Spe]1` files have **all** scaffolds, while the longer 
  `M_[species].p_ctg.hic_mapq0_noContigEC_minQ10_scaffolds_final_manualCuration.chr_M` has only one mitogenome (`chr_M`)
  with all scaffolds of alternative mitogenomes removed.

- `output/`
  - `minimap2/`: Output from `minimap2` rules.
    - `[CCS settings]/` eg `minPasses3_minRQ0.99/`: folder named after the settings used for CCS reads being mapped.
      
      Two classes of files here:
       - `[genome]-hifi.bam`: Reads from the genome indivdual (`mMyo[Spe]1`)to its consensus genome (`mMyo[Spe]1`).
       - `[ReferenceGenome]-[genome]-hifi.bam`: Reads from the genome individual mapped to an outgroup reference genome.

- `envs/`: contains all the environment YAML files for `conda`
- `logs/`: contains output logs
- `rules/`: contains all Snakemake rules, both test and final; finalized rules are imported in `./Snakemake`
- `pepsamples.tsv`: Portable Encapsulated Project (PEP) definition for HiFi samples
- `Snakefile`: master Snakefile that generates all the files in `output` 
