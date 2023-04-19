import os
import pandas as pd

def get_minimap_inputs(wildcards):
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    return ["data/HiFi-adapterFiltered/{species}/minPasses3_minRQ0.99/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq".format(species=s[0],
            pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]


rule minimap2_split:
    input:
       genome = "data/genomes/{species}/{contig}.fa",
       fastq = get_minimap_inputs
    #output: temp("output/minimap2-split/{species}/{contig}.unsorted.bam")
    output: "output/minimap2-split/{species}/{contig}.bam"
    conda: "../envs/minimap2.yaml"
    threads: 20
    shell: "minimap2 -t {threads} -ax map-hifi {input.genome} {input.fastq} | samtools sort -O BAM --reference {input.genome} -@ {threads} > {output}"
    #shell: "minimap2 -t {threads} -ax map-hifi {input.genome} {input.fastq} > {output}"

#rule minimap2_split_sort:
#    input:
#       genome = "data/genomes/{species}/{contig}.fa",
#       bam = "output/minimap2-split/{species}/{contig}.unsorted.bam"
#    output: "output/minimap2-split/{species}/{contig}.bam"
#    conda: "../envs/minimap2.yaml"
#    threads: 40
#    shell: "samtools sort -O BAM --reference {input.genome} -@ {threads} {input.bam} > {output}"


rule minimap2:
    input:
       genome = "data/genomes/{species}.fa",
       fastq = get_minimap_inputs
    output: "output/minimap2/{species}.bam"
    conda: "../envs/minimap2.yaml"
    threads: 32
    shell: "minimap2 -t {threads} -ax map-hifi {input.genome} {input.fastq} | samtools sort -O BAM --reference {input.genome} -@ {threads} > {output}"
