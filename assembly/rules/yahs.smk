rule samtools_sortByName:
    input: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.sorted.bam"
    output: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.sortedByReadName.bam"
    params:
        extra="-n",
    threads:  # Samtools takes additional threads through its option -@
        52     # This value - 1 will be sent to -@.
    wrapper:
        "v1.0.0/bio/samtools/sort"

rule yahs:
    input: 
        fasta = "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa",
        bam = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.sortedByReadName.bam"
    output: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_scaffolds_final.fa"
    params: 
        outdir = lambda wildcards: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}".format(species=wildcards.species, ccs_opts = wildcards.ccs_opts, hifiasm_opts = wildcards.hifiasm_opts, genometype = wildcards.genometype, hic = wildcards.hic)
    shell: "code/yahs/yahs {input.fasta} {input.bam} -o {params.outdir}"
