Assembly of Chromosome-Level Genomes
====================================


Methods:
-------

Our genomes were assembled using PacBio CCS reads, with Dovetail Omni-C sequencing for phasing 
and scaffolding. Briefly, primary cell lines were derived from 3mm wing punches of wild-caught 
bats, and expanded in culture to ~20M cells per line. Cell pellets were then collected for either 
DNA extraction or for Hi-C library prep. The PacBio reads were processed using SMRTTools (version 
TODO, @smrttools) to generate the circular consensus sequences using the settings 
`--minPasses=3 --minRQ=0.99`. Hi-C reads were processed using `trimmomatic` (version TODO, 
@trimmomatic) to remove adapter sequences and low-quality bases using the settings 
`ILLUMINACLIP:data/trimmomatic-adapters/TruSeq3-PE-2.fa:2:40:15 SLIDINGWINDOW:5:20`.
To generate the primary contig assemblies, we used `hifiasm` (version TODO, @hifiasmhic) in Hi-C 
mode, providing both the CCS reads and the trimmed Hi-C reads as input, and removing all haplotigs 
with `-l2`. For our reference genomes, we proceeded with the primary contig assembly 
(`*.asm.hic.p_ctg.gfa`).

To scaffold our reference genomes, we used YAHS (version TODO, @yahs) with the primary contig 
assembly and our Hi-C data. Dovetail Omni-C data were processed and mapped to the genome following 
the manufacturer's instructions using `bwa` (version TODO, @bwa), `pairtools` (version TODO, 
@pairtools), and `samtools` (version TODO, @samtools). YAHS was run using both default settings as
well as with `--no-contig-ec`; after comparing the outputs, we proceeded with the `--no-contig-ec` 
version for our final assemblies. 

To finalize the assemblies, we performed manual curation using PreTextView and the Rapid Curation 
toolkit (version TODO, @pretext ; @rapidcuration). Sex chromosomes were identified based on half-
coverage in XY genomes, and the X chromsomes were further identified after comparisons to XX 
genomes. Mitochondrial genomes were identified and removed from the final assembly by running 
`mitohifi` (version 2.TODO, @mitohifi) in contig mode on the assembly and removing all scaffolds
identified as mitogenomes. The consensus mitogenome from `mitohifi` was designated as the 
representative mitogenome for the assembly. 


Repo Structure: 
---------------

- `data/`: Contains all initial data for assembly, including sequencing data. 
Also contains folders with symlinks to final products;
- `envs/`: Contains the Conda environment files defining the tools for assembly. 
Env files with dates are the explicit environments used in generating the assemblies;
- `output/`: Contains all the intermediate and final outputs of assembly;
- `rules/`: Contains the Snakefiles with rules for assembling and QC-ing the genomes;
- `Snakefile`: The master Snakefile for assembling the genomes;
- `pepsamples.tsv`: Table with all the sequencing metadata (derived from the files themselves);
- `pep_mitochondrion.tsv`: Table listing representative mitogenomes from NCBI used with `mitohifi`.
- `rna_pepsamples.tsv`: Table with all the sequencing metadata for the RNA-seq samples in this study.
