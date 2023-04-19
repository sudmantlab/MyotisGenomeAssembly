rule get_flip_table:
    version: 0.1
    input: "output/wfmash/1-to-1/{ref_genome}_{query_genome}.approx.paf"
    output: "output/flip_table/{ref_genome}_{query_genome}.csv"
    conda: "../envs/r-pafr.yaml"
    script: "../code/get_flip_table.R"


rule flip_genome:
    version: 0.1
    input: 
        genome = "data/genomes/{genome}.fa",
        flip_table = "output/flip_table/{ref_genome}_{genome}.csv"
    output: "output/flipped_genomes/{genome}_rel-{ref_genome}.fa"
    threads: 40
    conda: "../envs/biopython.yaml"
    script: "../code/flip_chr.py"

rule rename_genome:
    version: 0.1
    input: 
        genome = "data/genomes/{genome}.fa",
        flip_table = ancient("output/flip_table/{ref_genome}_{genome}.csv")
    output: "output/renamed_flipped_genomes/{genome}_rel-{ref_genome}_pangenomeNamed.fa"
    threads: 40
    conda: "../envs/biopython.yaml"
    script: "../code/rename_chr.py"

