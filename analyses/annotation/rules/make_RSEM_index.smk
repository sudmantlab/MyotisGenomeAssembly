# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import tempfile
import pandas as pd
import csv
import os
import pdb

#def get_inputs(wildcards):
#    inputs = []
#    for species, ref in config['reference_by_species'].items():
#        inputs.append("output/RSEM_indexes/"
#                      "{annot}/"
#                      "{species}/"
#                      "{ref}/"
#                      "{ref}.ti".format(annot=config['annotation'],
#                                        species=species,
#                                        ref=ref))
#        inputs.append("output/RSEM_indexes/"
#                      "{annot}/"
#                      "{species}/"
#                      "{ref}/"
#                      "{ref}.ngvec".format(annot=config['annotation'],
#                                           species=species,
#                                           ref=ref))
#    return inputs

#rule rsem_index_all:
#    input:
#        get_inputs

rule rsem_make_ngvec:
    output:
        #"output/RSEM_indexes/funannotate/{species}/{ref}/{ref}.ngvec"
        "output/RSEM_indexes/{species}/{ref}-{type}-{id}/{ref}.ngvec"
    input:
        ti = "output/RSEM_indexes/{species}/{ref}-{type}-{id}/{ref}.ti",
        transcripts = "output/RSEM_indexes/{species}/{ref}-{type}-{id}/{ref}.transcripts.fa"
    threads: 1
    conda: '../envs/STAR-RSEM-EBSeq.yaml'
    shell: "perl code/rsem-generate-ngvector "
           "{input.transcripts} "
           "output/RSEM_indexes/{wildcards.species}/{wildcards.ref}-{wildcards.type}-{wildcards.id}/{wildcards.ref}"


rule rsem_prepare_reference: 
    input:
        #gff = 'output/funannotate/TOGA_Prot/{ref}_TOGA_Prot/predict_results/{species}.gff3',
        gff = "data/GFF_evidences/{type}/{ref}_{id}.gff3",
        fasta = 'data/genomes/{ref}.fa'
    output: 
        "output/RSEM_indexes/{species}/{ref}-{type}-{id}/{ref}.ti",
        "output/RSEM_indexes/{species}/{ref}-{type}-{id}/{ref}.transcripts.fa"
    threads: 1
    conda: '../envs/STAR-RSEM-EBSeq.yaml'
    shell: "rsem-prepare-reference "
           " --gff3 {input.gff} "
           " {input.fasta} "
           " output/RSEM_indexes/{wildcards.species}/{wildcards.ref}-{wildcards.type}-{wildcards.id}/{wildcards.ref} "

