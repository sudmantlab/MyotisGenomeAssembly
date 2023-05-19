rule samtools_faidx_yahs:
    input: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa"
    output: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa.fai"
    conda: "../envs/omni-c.yaml"
    shell: "samtools faidx {input}"
    

rule genome_file_yahs:
    input: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa.fai"
    output: "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.genome"
    shell: "cut -f1,2 {input} > {output}"


rule bwa_index_yahs:
    input:
        "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa",
    output:
        idx=multiext("output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{hic}.log",
    conda: "../envs/omni-c.yaml"
    shell: "bwa index {input}"


rule bwa_aln_yahs:
    input:
        reads= get_omnic_reads,
        # Index can be a list of (all) files created by bwa, or one of them
        fasta="output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa",
        idx=multiext("output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.aligned.sam",
    params:
        extra="-5SP -T0"
    log:
        "logs/bwa_aln/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.aligned.log",
    threads: 32
    conda: "../envs/omni-c.yaml"
    shell: "bwa mem -5SP -T0 -t{threads} {input.fasta} {input.reads} -o {output}"
    #wrapper:
    #    "v1.1.0/bio/bwa/mem"


rule pairtools_parse_yahs:
    input: 
        sam = "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.aligned.sam",
        genome = "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.genome"
    output: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.aligned.pairsam"
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.aligned.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: "pairtools parse "
           "--min-mapq 40 "
           "--walks-policy 5unique "
           "--max-inter-align-gap 30 "
           "--nproc-in {threads} "
           "--nproc-out {threads} "
           "--chroms-path {input.genome} "
           "{input.sam} > {output} 2> {log}"


#rule pairtools_sort_yahs:
#    input: 
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.aligned.pairsam"
#    output: 
#        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.sorted.pairsam")
#    log: 
#        "logs/pairtools_sort/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.sorted.log"
#    conda: "../envs/omni-c.yaml"
#    threads: 32
#    shell: 
#        """
#        mkdir /global/scratch2/mvazquez/tmp_pairtools
#        pairtools sort \
#          --tmpdir=/global/scratch2/mvazquez/tmp_pairtools \
#          --nproc {threads} \
#          {input} > {output} 2> {log}
#        rm -r /global/scratch2/mvazquez/tmp_pairtools
#        """


#rule pairtools_dedup_yahs:
#    input: 
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.sorted.pairsam"
#    output: 
#        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.dedup.pairsam")
#    log: 
#        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.dedup.log"
#    conda: "../envs/omni-c.yaml"
#    threads: 32
#    shell: 
#        """
#        pairtools dedup \
#         --nproc-in {threads} \
#         --nproc-out {threads} \
#         --mark-dups \
#         --output-stats {output}.stats \
#         --output {output} {input} 2> {log}
#        """


#rule pairtools_split_yahs:
#    input: 
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.dedup.pairsam"
#    output: 
#        bam = temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.bam"),
#        pairs = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.pairs",
#    log: 
#        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.split.log"
#    conda: "../envs/omni-c.yaml"
#    threads: 32
#    shell: 
#        """
#        pairtools split \
#        --nproc-in {threads} \
#        --nproc-out {threads} \
#        --output-pairs {output.pairs} \
#        --output-sam {output.bam} \
#        {input} 2> {log}
#        """


#rule samtools_sort_yahs:
#    input:
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.bam"
#    output:
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.sorted.bam"
#    params:
#        #extra = "-m 100G"
#    threads:  # Samtools takes additional threads through its option -@
#        32     # This value - 1 will be sent to -@.
#    wrapper:
#        "v1.0.0/bio/samtools/sort"


#rule samtools_index_yahs:
#    input:
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.sorted.bam"
#    output:
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.sorted.bam.bai"
#    log:
#        "output/Omni-C_samtools_index/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}_scaffolds_final.sorted.log"
#    params:
#        "" # optional params string
#    threads:  # Samtools takes additional threads through its option -@
#        32     # This value - 1 will be sent to -@
#    wrapper:
#        "0.77.0/bio/samtools/index"
