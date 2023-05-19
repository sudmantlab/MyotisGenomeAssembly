import os
import tempfile
import pandas as pd
import pdb

#configfile: "config.yaml"


def get_STAR_all(wildcards):
    inputs = []
    out_path= "output/mapping/{species}/{genome}/{sample}/Aligned.sortedByCoord.out.bam.bai"
    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
    for species in config['species']:
        #print(species)
        for ref in config['reference_by_species'][species]:
            #print(ref)
            inputs.extend([out_path.format(species = species, genome = ref, sample = r[3]) for r in samples[samples["Species"] == species].itertuples()])
            #print(inputs)
    return inputs


rule STAR_all:
    input:
        get_STAR_all
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem-per-cpu 1000 "
                                     "--mem 4000 "
                                     "--time 0-0:05:00 "
                                     "-J rule_all "
                                     "-o logs/rule_all_%j.logs "
                                     "-p defq "

rule index_bams:
    input:
        "output/mapping/{species}/{genome}/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "output/mapping/{species}/{genome}/{sample}/Aligned.sortedByCoord.out.bam.bai"
    threads: 13
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 30000 "
                                     "--time 0-6:00:00 "
                                     "-J index_{sample} "
                                     "-o logs/index_{sample}_%j.logs "
                                     "-p defq "
                                     "".format(sample=wildcards.sample)
    shell:
        "samtools index -@ {threads} {input[0]}"


def get_map_fastq_inputs(wildcards):
    path_trimmed = "data/RNA-seq_trimmed/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    import pandas as pd
    from pathlib import Path
    samples = pd.read_table("rna_pepsamples.tsv", index_col= False)
    samples = samples[samples.Sample == wildcards.sample]
    input_dict = {"left": [], "right": [], "index": "output/STAR_indexes/{species}/{genome}/SA"}
    samples_grouped = samples.groupby(samples.Sample)
    for sample in set(samples["Sample"].tolist()):
        sample_subset = samples_grouped.get_group(sample)
        #print(sample_subset)
        if len(sample_subset) == 0:
            raise Exception("No files available for sample {}".format(sample))
        #print(r[2])
        input_dict["left"].extend([path_trimmed.format(species = config['species_abbr'][r[2]], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["Read"] == "R1"].itertuples()])
        input_dict["right"].extend([path_trimmed.format(species = config['species_abbr'][r[2]], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["Read"] == "R2"].itertuples()])
    print(input_dict)
    return input_dict


rule map_pe:
    input: unpack(get_map_fastq_inputs)
    output:
        "output/mapping/{species}/{genome}/{sample}/Aligned.sortedByCoord.out.bam",
        "output/mapping/{species}/{genome}/{sample}/Log.final.out"
    shadow: "shallow"
    log: "logs/mapping/{species}/{genome}/{sample}/star.log"
#    params:
#        map_opts = config['map_opts']
    conda: "../envs/STAR.yaml"
    threads: 40
    shell:
           "STAR runThreadN {threads} "
           "--genomeDir output/STAR_indexes/{wildcards.species}/{wildcards.genome} "
           "--readFilesIn {input.left} {input.right} "
           "--readFilesCommand zcat "
           "{params.map_opts} "
           "--outFileNamePrefix output/mapping/{wildcards.species}/{wildcards.genome}/{wildcards.sample}/ "
           "--outSAMtype BAM SortedByCoordinate  "
           "--twopassMode Basic "
           "--limitBAMsortRAM 60000000000 "
           "> {log} 2>&1"

rule clean:
    shell:
        "rm -rf output/mapping; rm logs/*; rm -rf output/rsem"

