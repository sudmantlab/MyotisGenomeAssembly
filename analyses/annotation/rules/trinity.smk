localrules: trinityStats

rule trinity_pe:
    input: 
        left="output/trimmed/{species}/PAIRED/{sample}-trimmed_R1.fastq.gz",
        right="output/trimmed/{species}/PAIRED/{sample}-trimmed_R2.fastq.gz"
    output:
        "output/trinity/{species}/PAIRED/{sample}-trinity/Trinity.fasta"
    log:
        'logs/trinity/{species}/PAIRED/{sample}-trinity.log'
    params:
        extra="--SS_lib_type=RF"
    threads: 40
    resources:
        mem_gb=192
    wrapper:
        "v1.14.0/bio/trinity"

rule trinity_SE:
    input: 
        left="output/trimmed/{species}/SINGLE/{sample}-trimmed.fastq.gz"
    output:
        "output/trinity/{species}/SINGLE/{sample}-trinity/Trinity.fasta"
    log:
        'logs/trinity/{species}/SINGLE/{sample}-trinity.log'
    threads: 40
    resources:
        mem_gb=192
    wrapper:
        "v1.14.0/bio/trinity"


rule trinityStats: 
    input: "output/trinity/{species}/{paired}/{sample}-trinity/Trinity.fasta"
    output: "output/trinity/{species}/{paired}/{sample}-trinity/TrinityStats.txt"
    shell: "TrinityStats.pl {input} > {output}"

