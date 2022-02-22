rule samtools_sortname_fixmate_sortpos_rmdup:
    input: 'output/mapping/{sample}/Aligned.sortedByCoord.out.bam'
    output:
        bam = temp('output/bam_rmdup/bam/{sample}.noDup.bam'),
        stats = 'output/bam_rmdup/stats/{sample}.txt'
    conda: "../envs/STAR-EBSeq-RSEM.yaml"
    threads: 20
    #benchmark:
    #    repeat("benchmarks/samtools-piped/{sample}.tsv",3)
    shadow: "shallow"
    shell: """
      mkdir -p output/bam_rmdup/{{bam,stats}};
      samtools sort -n -@ {threads} -O BAM -o - {input} | 
        samtools fixmate -m -@ {threads} -O BAM - - | 
        samtools sort -@ {threads} -O BAM -o - - | 
        samtools markdup -rs -f {output.stats} -O BAM -@ {threads} - {output.bam}
    """

# WASP

rule WASP_samtools_sortname_fixmate_sortpos_rmdup:
    input: 'output/mapping-WASP/{sample}/Aligned.sortedByCoord.out.bam'
    output:
        bam = temp('output/WASP_bam_rmdup/bam/{sample}.noDup.bam'),
        stats = 'output/bam_rmdup/stats/{sample}.txt'
    conda: "../envs/STAR-EBSeq-RSEM.yaml"
    threads: 20
    #benchmark:
    #    repeat("benchmarks/samtools-piped/{sample}.tsv",3)
    shadow: "shallow"
    shell: """
      mkdir -p output/bam_rmdup/{{bam,stats}};
      samtools sort -n -@ {threads} -O BAM -o - {input} | 
        samtools fixmate -m -@ {threads} -O BAM - - | 
        samtools sort -@ {threads} -O BAM -o - - | 
        samtools markdup -rs -f {output.stats} -O BAM -@ {threads} - {output.bam}
    """


