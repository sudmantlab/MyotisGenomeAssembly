rule samtools_faidx:
    input: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.fa"
    output: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.fa.fai"
    conda: "../envs/omni-c.yaml"
    shell: "samtools faidx {input}"
    

rule genome_file:
    input: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.fa.fai"
    output: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.genome"
    shell: "cut -f1,2 {input} > {output}"


rule bwa_index:
    input:
        "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.fa",
    output:
        idx=multiext("output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{hic}.log",
    params:
        #algorithm="bwtsw",
    conda: "../envs/omni-c.yaml"
    shell: "bwa index {input}"
    #wrapper:
    #    "v1.2.0/bio/bwa/index"

def get_omnic_reads(wildcards):
    path_trimmed = "output/trimmed-hic/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    input_list = []
    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
    samples = samples[samples.species==wildcards.species][samples.type == "Hi-C"]
    samples_grouped = samples.groupby(samples.sample_name)
    for sample in set(samples["sample_name"].tolist()):
        sample_subset = samples_grouped.get_group(sample)
        if len(sample_subset) == 0:
            raise Exception("No files available for sample {}".format(sample))
        input_list.extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R1"].itertuples()])
        input_list.extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R2"].itertuples()])
    return input_list


rule bwa_aln:
    input:
        reads= get_omnic_reads,
        # Index can be a list of (all) files created by bwa, or one of them
        fasta="output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.fa",
        idx=multiext("output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.aligned.sam",
    params:
        extra="-5SP -T0"
    log:
        "logs/bwa_aln/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.aligned.log",
    threads: 32
    conda: "../envs/omni-c.yaml"
    shell: "bwa mem -5SP -T0 -t{threads} {input.fasta} {input.reads} -o {output}"
    #wrapper:
    #    "v1.1.0/bio/bwa/mem"


rule pairtools_parse:
    input: 
        sam = "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.aligned.sam",
        genome = "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.genome"
    output: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.aligned.pairsam"
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.aligned.log"
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


rule pairtools_sort:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.aligned.pairsam"
    output: 
        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.pairsam")
    log: 
        "logs/pairtools_sort/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        mkdir /global/scratch2/mvazquez/tmp_pairtools
        pairtools sort \
          --tmpdir=/global/scratch2/mvazquez/tmp_pairtools \
          --nproc {threads} \
          {input} > {output} 2> {log}
        rm -r /global/scratch2/mvazquez/tmp_pairtools
        """


rule pairtools_dedup:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.pairsam"
    output: 
        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.dedup.pairsam")
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.dedup.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        pairtools dedup \
         --nproc-in {threads} \
         --nproc-out {threads} \
         --mark-dups \
         --output-stats {output}.stats \
         --output {output} {input} 2> {log}
        """


rule pairtools_split:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.dedup.pairsam"
    output: 
        bam = temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.bam"),
        pairs = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs",
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.split.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        pairtools split \
        --nproc-in {threads} \
        --nproc-out {threads} \
        --output-pairs {output.pairs} \
        --output-sam {output.bam} \
        {input} 2> {log}
        """


rule samtools_sort:
    input:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.bam"
    output:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.bam"
    params:
        #extra = "-m 100G"
    threads:  # Samtools takes additional threads through its option -@
        32     # This value - 1 will be sent to -@.
    wrapper:
        "v1.0.0/bio/samtools/sort"


rule samtools_index:
    input:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.bam"
    output:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.bam.bai"
    log:
        "output/Omni-C_samtools_index/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        32     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"
