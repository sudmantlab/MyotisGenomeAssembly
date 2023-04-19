# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Make a STAR index
Expects a json config file with the following structure, 
   

{
    "reference_by_species": {"homo_sapiens": "hg19",
                           "bos_taurus" : "bosTau6",
                           "macaca_mulatta" : "rheMac2",
                           "mus_musculus" : "mm10",
                           "gallus_gallus" : "galGal4",
                           "rattus_norvegicus" : "rn5"},
    "genome_path": "/scratch/users/psudmant/genomes/UCSC/UCSC/"
}
"""

__author__ = "Peter Sudmant"
__license__ = "MIT"

#configfile: "config_trimmed.json"

import tempfile

def get_inputs(wildcards):
    inputs = []
    for species, ref in config['reference_by_species'].items():
        inputs.append("STAR_indexes/{species}/{ref}/SA".format(species=species,
                                                               ref=ref))
    return inputs

#rule all:
#    input:
#        get_inputs
#    params:
#        slurm_opts=lambda wildcards: "-n1 "
#                                     "--share "
#                                     "--export ALL "
#                                     "--mem-per-cpu 1000 "
#                                     "--mem 4000 "
#                                     "--time 0-8:00:00 "
#                                     "-J rule_all "
#                                     "-p defq "
#                                     "-o logs/rule_all_%j.logs "


rule STAR_index_creation:
    input:
        lambda wildcards: "{genome_path}/"
                          "{ref}.fa".format(genome_path=config['genome_path'],
                                            #species=wildcards.species,
                                            ref=config['reference_by_species'][wildcards.species])
    output:
        "output/STAR_indexes/{species}/{ref}/SA"
    threads: 
        64 
    # conda: "envs/STAR.yml"
    params:
        slurm_opts=lambda wildcards: "-n24 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 100000 "
                                     "-m cyclic "
                                     "--time 0-8:00:00 "
                                     "-o logs/STAR_idx_{species}_%j.logs "
                                     "-J STAR_idx_{species} "
                                     "-p defq ".format(species=wildcards.species)
    run:
        tmp_dir = tempfile.TemporaryDirectory(prefix=wildcards.species,dir="/dev/shm")
        
        shell("STAR --runThreadN {threads} "
              "--runMode genomeGenerate "
              "--genomeDir {tmp_dir} "
              "--genomeFastaFiles {input_fa} "
              "--outFileNamePrefix logs/STAR_index_{species} "
              "--limitGenomeGenerateRAM=124544990592 "
              "".format(input_fa=input[0],
                        tmp_dir=tmp_dir.name,
                        threads=threads,
                        species=wildcards.species))
        shell("cp -r {tmp_path}/* output/STAR_indexes/{species}/{ref} "
                             " ".format(tmp_path=tmp_dir.name,
                                        species=wildcards.species,
                                        ref = wildcards.ref))
        """cleanup the tempdir"""
        tmp_dir.cleanup()

