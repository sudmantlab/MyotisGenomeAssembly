import os
import pandas as pd 

def get_freebayes_inputs(wildcards):
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    return ["output/minimap2/{species}/{pacbio1}/{pacbio2}/{id}.sam".format(species=s[0],
            pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]

rule samtools_index_split:
    input:
        "output/minimap2-split/{species}/{contig}.bam"
    output:
        "output/minimap2-split/{species}/{contig}.bam.bai"
    log:
        "logs/minimap2-split/{species}/{contig}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"

rule freebayes_contig_split:
    input:
        ref="data/genomes/{species}/{contig}.fa",
        # you can have a list of samples here
        samples= "output/minimap2-split/{species}/{contig}.bam",
        # the matching BAI indexes have to present for freebayes
        indexes= "output/minimap2-split/{species}/{contig}.bam.bai",
    output:
        "output/freebayes-calls-split/{species}-{contig}.vcf"  # either .vcf or .bcf
    log:
        "logs/freebayes/{species}-{contig}.log"
    params:
        extra="-! 5 -F 0.5",         # optional parameters
        # insert sisze for pacbio is
        #          15000
        chunksize=150000, # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # flag to use bcftools norm to normalize indels
    threads: 32
    wrapper:
        "v0.75.0/bio/freebayes"


rule getHetSites:
    input: "output/freebayes-calls-split/{species}-{contig}.vcf"
    output: "output/heterozygousSites-split/{species}/{species}-{contig}.bed"
    threads: 1
    conda: "../envs/pyvcf.yaml"
    shell: "python code/getHet.py {input} {output}"


rule getHetSites_joint:
    input: "output/freebayes-calls/{species}.vcf"
    output: "output/heterozygousSites/{species}.bed"
    threads: 1
    conda: "../envs/pyvcf.yaml"
    shell: "python code/getHet.py {input} {output}"


def get_contigs(wildcards):
    with open("data/genomes/{species}.fa.fai".format(species=wildcards.species)) as infile:
        return [line.split("\t")[0] for line in infile]

def get_merge_het_sites(wildcards):
    return expand("output/heterozygousSites-split/{species}/{species}-{contig}.bed", contig = get_contigs(wildcards), species = wildcards.species)

rule mergeHetSites:
    input: get_merge_het_sites
    output: "output/heterozygousSites/{species}.bed"
    shell: "cat {input} > {output}"

rule bedtools_makewindows:
    input: "data/genomes/{species}.fa.fai",
    output: "data/BED_windows/{species}.{nKB}KBWindow_{nKB2}KBSteps.bed"
    conda: "../envs/bedtools.yaml"
    shell: "bedtools makewindows -g {input} -w {wildcards.nKB}000 -s {wildcards.nKB2}000 > {output}"

rule bedtools_coverage:
    input: 
        a = "data/BED_windows/{species}.{nKB}KBWindow_{nKB2}KBSteps.bed",
        b = "output/heterozygousSites/{species}.bed"
    output: "output/heterozygosity_windowed/{species}/{species}.{nKB}KBWindow_{nKB2}KBSteps.bed"
    conda: "../envs/bedtools.yaml"
    shell: "bedtools coverage -a {input.a} -b {input.b} > {output}"
    

rule bgzip_vcf:
    input: "output/freebayes-calls-split/{species}-{contig}.vcf"
    output: 
        gz="output/freebayes-calls-split/{species}-{contig}.vcf.gz",
        gzi="output/freebayes-calls-split/{species}-{contig}.vcf.gz.gzi"
    threads: 32
    shell: "bgzip -@ {threads} -if {input}"
