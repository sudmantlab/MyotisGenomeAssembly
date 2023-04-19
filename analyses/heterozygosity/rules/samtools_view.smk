rule samtools_index:
    input:
        "output/minimap2/{species}.bam"
    output:
        "output/minimap2/{species}.bam.bai"
    log:
        "logs/minimap2/{species}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"

rule samtools_view:
    input: "output/minimap2/{species}.bam"
    output: "output/minimap2-split/{species}/{contig}.bam"
    shell: "samtools view -b {input} {wildcards.contig} > {output}"

