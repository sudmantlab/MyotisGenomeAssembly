
def get_seqinfo(wildcards):
    import pandas as pd
    samples = pd.read_table(config["input_seqinfo"], index_col="fastq")
    #print(samples)
    #print(wildcards.fastq)
    samples = samples.filter(items=[wildcards.fastq], axis=0)
    samples = samples.to_dict('records')[0]
    #print(samples)
    samples["sample_id"] = config["sample_experiment"][wildcards.fastq]
    readgroup = '-r ID:{s_id}_{f_id}_{f_lane} -r SM:{s_id} -r LB:{lib} -r PU:{f_lane} -r PL:ILLUMINA'.format(f_id=samples["flowcell_id"],
                                                                                                             f_lane=samples["flowcell_lane"],
                                                                                                             lib=samples["library"],
                                                                                                             s_id=samples["sample_id"])
    #print(readgroup)
    return readgroup

rule add_readgroup:
    input: 'output/bam_rmdup/bam/{fastq}.noDup.bam'
    output: 'output/bam_rmdup/bam/{fastq}.noDup.rg.bam'
    params: 
        readgroup = get_seqinfo
    conda: "../envs/STAR-EBSeq-RSEM.yaml"
    shell:
        'samtools addreplacerg {params.readgroup} -o {output} {input}'


rule index_bams_rg:
    input:
        "output/bam_rmdup/bam/{fastq}.noDup.rg.bam"
    output:
        "output/bam_rmdup/bam/{fastq}.noDup.rg.bam.bai"
    threads: 40
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 30000 "
                                     "--time 0-6:00:00 "
                                     "-J index_{sample} "
                                     "-o logs/index_{sample}_%j.logs "
                                     "-p defq "
                                     "".format(sample=wildcards.fastq)
    shell:
        "samtools index -@ {threads} {input}"


# WASP
rule add_readgroup_WASP:
    input: 'output/WASP_bam_rmdup/bam/{fastq}.noDup.bam'
    output: 'output/WASP_bam_rmdup/bam/{fastq}.noDup.rg.bam'
    params:
        readgroup = get_seqinfo
    conda: "../envs/STAR-EBSeq-RSEM.yaml"
    shell:
        'samtools addreplacerg {params.readgroup} -o {output} {input}'


rule index_bams_rg_WASP:
    input:
        "output/WASP_bam_rmdup/bam/{fastq}.noDup.rg.bam"
    output:
        "output/WASP_bam_rmdup/bam/{fastq}.noDup.rg.bam.bai"
    threads: 40
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 30000 "
                                     "--time 0-6:00:00 "
                                     "-J index_{sample} "
                                     "-o logs/index_{sample}_%j.logs "
                                     "-p defq "
                                     "".format(sample=wildcards.fastq)
    shell:
        "samtools index -@ {threads} {input}"
