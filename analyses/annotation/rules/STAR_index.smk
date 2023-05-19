__author__ = "Peter Sudmant"
__license__ = "MIT"

#configfile: "config.yaml"

import tempfile

#def get_inputs(wildcards):
#    inputs = []
#    for species, ref in config['reference_by_species'].items():
#        inputs.extend(["output/STAR_indexes/{species}/{ref}/SA".format(species=species, ref=r) for r in ref])
#    return inputs

#rule all:
#    input:
#        get_inputs

rule STAR_index_creation:
    input:
        fasta = "data/genomes/{ref}.fa",
        #gff = "output/funannotate/TOGA_Prot/mMyoOcc1_TOGA_Prot/predict_results/{species}.gff"
    output:
        multiext("output/STAR_indexes_noGFF/{species}/{ref}/", "chrLength.txt", "chrNameLength.txt", "chrName.txt", "chrStart.txt", "Genome", "genomeParameters.txt", "SA", "SAindex")
    threads: 
        32 
    shadow: "shallow"
    conda: "../envs/STAR-RSEM-EBSeq.yaml"
    shell: "STAR --runThreadN {threads} "
           "--runMode genomeGenerate "
           #"--sjdbGTFtagExonParentTranscript Parent "
           #"--sjdbGTFfile {input.gff} "
           "--genomeDir output/STAR_indexes/{wildcards.species}/{wildcards.ref} "
           "--genomeFastaFiles {input.fasta} "
           "--outFileNamePrefix logs/STAR_index/{wildcards.species}-{wildcards.ref} "
           "--limitGenomeGenerateRAM=124544990592 "

rule STAR_index_creation_withGFF:
    input:
        fasta = "data/genomes/{ref}.fa",
        gff = "output/funannotate/TOGA_Prot/{ref}_TOGA_Prot/predict_results/{species}.gff3"
    output:
        multiext("output/STAR_indexes/{species}/{ref}/", "chrLength.txt", "chrNameLength.txt", "chrName.txt", "chrStart.txt", "Genome", "genomeParameters.txt", "SA", "SAindex")
    threads: 
        32 
    shadow: "shallow"
    conda: "../envs/STAR-RSEM-EBSeq.yaml"
    shell: "STAR --runThreadN {threads} "
           "--runMode genomeGenerate "
           "--sjdbGTFtagExonParentTranscript Parent "
           "--sjdbGTFfile {input.gff} "
           "--genomeDir output/STAR_indexes/{wildcards.species}/{wildcards.ref} "
           "--genomeFastaFiles {input.fasta} "
           "--outFileNamePrefix logs/STAR_index/{wildcards.species}-{wildcards.ref} "
           "--limitGenomeGenerateRAM=124544990592 "

