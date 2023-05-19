import os
import pandas as pd 

localrules: collate_bams

def get_all_minimap2_species(wildcards):
    f_path = "output/minimap2/{ccs_settings}/{species}/{pacbio1}/{pacbio2}/{id}.sorted.bam"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    input_samples = [f_path.format(species=s[0], ccs_settings = wildcards.ccs_settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    print(input_samples[0])
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples

def get_all_minimap2_species_distinctRef(wildcards):
    f_path = "output/minimap2/{ccs_settings}/{Refspecies}/{species}/{pacbio1}/{pacbio2}/{id}.sorted.bam"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    input_samples = [f_path.format(species=s[0], Refspecies = wildcards.Refspecies, ccs_settings = wildcards.ccs_settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    print(input_samples[0])
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples

def get_seqinfo_wc(wildcards):
    machine_date_time_tuple = wildcards.machine_date_time.split("_")
    rg_dict = dict(
                   f_id = wildcards.machine_date_time + "_" + wildcards.cell + "_" + wildcards.id,
                   s_id = wildcards.id,
                   f_lane = wildcards.cell
                   )
    readgroup = '@RG\\tID:{f_id}\\tSM:{s_id}\\tPU:{f_lane}\\tPL:PACBIO'.format(**rg_dict)
    print(readgroup)
    return readgroup

rule minimap2_bam: 
    version: "0.3"
    input:
        refgenome = "data/genomes/{species}.fa",
        hifi = "data/HiFi-adapterFiltered/{species}/{ccs_settings}/{machine_date_time}/{cell}/{id}.ccs.filt.fastq.gz"
    output: temp("output/minimap2/{ccs_settings}/{species}/{machine_date_time}/{cell}/{id}.bam")
    params:
        readgroup= get_seqinfo_wc
    conda: "../envs/minimap2.yaml"
    threads: 40
    shell: "minimap2 --version && minimap2 -R '{params.readgroup}' -t {threads} -ax map-hifi {input.refgenome} {input.hifi} | samtools view -b > {output}"

rule minimap2_bam_distinctRef: 
    version: "0.3"
    input:
        refgenome = "data/genomes/{Refspecies}.fa",
        hifi = "data/HiFi-adapterFiltered/{species}/{ccs_settings}/{machine_date_time}/{cell}/{id}.ccs.filt.fastq.gz"
    output: temp("output/minimap2/{ccs_settings}/{Refspecies}/{species}/{machine_date_time}/{cell}/{id}.bam")
    params:
        readgroup= get_seqinfo_wc
    conda: "../envs/minimap2.yaml"
    threads: 40
    shell: "minimap2 --version && minimap2 -R '{params.readgroup}' -t {threads} -ax map-hifi {input.refgenome} {input.hifi} | samtools view -b > {output}"


rule samtools_sort:
    version: "0.1"
    input: "output/minimap2/{ccs_settings}/{species}/{machine_date_time}/{cell}/{sample}.bam"
    output: temp("output/minimap2/{ccs_settings}/{species}/{machine_date_time}/{cell}/{sample}.sorted.bam")
    threads: 32
    conda: "../envs/minimap2.yaml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

#rule samtools_sort_distinctRef:
#    version: "0.1"
#    input: "output/minimap2/{ccs_settings}/{Refspecies}/{species}/{machine_date_time}/{cell}/{id}.bam"
#    output: temp("output/minimap2/{ccs_settings}/{Refspecies}/{species}/{machine_date_time}/{cell}/{id}.sorted.bam")
#    threads: 32
#    conda: "../envs/minimap2.yaml"
#    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    version: "0.1"
    input: get_all_minimap2_species
    output: "output/minimap2/{ccs_settings}/{species}-hifi.bam"
    threads: 32
    conda: "../envs/minimap2.yaml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input} "

rule collate_bams_distinctRef:
    version: "0.1"
    input: get_all_minimap2_species_distinctRef
    output: "output/minimap2/{ccs_settings}/{Refspecies}-{species}-hifi.bam"
    threads: 32
    conda: "../envs/minimap2.yaml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input} "
