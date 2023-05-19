rule samtools_sortByName:
    input: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam"
    output: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sortedByReadName.bam"
    params:
        extra="-n",
    threads:  # Samtools takes additional threads through its option -@
        52     # This value - 1 will be sent to -@.
    wrapper:
        "v1.0.0/bio/samtools/sort"

rule yahs:
    input: 
        fasta = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa",
        bam = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sortedByReadName.bam"
    output: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}_scaffolds_final.fa"
    params: 
        outdir = "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}"
    shell: "code/yahs/yahs {input.fasta} {input.bam} -o {params.outdir}"


rule yahs_noContigEC:
    input: 
        fasta = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa",
        bam = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sortedByReadName.bam"
    output: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}_noContigEC_minQ{minmapq}_scaffolds_final.fa"
    params: 
        outdir = "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}_noContigEC_minQ{minmapq}"
    shell: "code/yahs/yahs {input.fasta} {input.bam} -o {params.outdir} --no-contig-ec -q {wildcards.minmapq}"


rule link_yahs_assembly:
    input: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}_scaffolds_final.fa"
    output: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}_scaffolds_final.fa"
    shell: "ln -s $(pwd -P)/{input} {output}"

rule link_yahs_EC_assembly:
    input: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}_noContigEC_minQ{minmapq}_scaffolds_final.fa"
    output: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_mapq{mapq}_noContigEC_minQ{minmapq}_scaffolds_final.fa"
    shell: "ln -s $(pwd -P)/{input} {output}"
 
