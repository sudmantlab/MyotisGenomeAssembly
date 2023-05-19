rule samtools_faidx:
    input: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa"
    output: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa.fai"
    conda: "../envs/omni-c.yaml"
    shell: "samtools faidx {input}"
    

rule genome_file:
    input: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa.fai"
    output: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.genome"
    shell: "cut -f1,2 {input} > {output}"


rule bwa_index:
    input:
        "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa",
    output:
        idx=multiext("data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.log",
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
        fasta="data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa",
        idx=multiext("data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        temp("output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.aligned.sam")
    params:
        extra="-5SP -T0"
    log:
        "logs/bwa_aln/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.aligned.log",
    threads: 32
    conda: "../envs/omni-c.yaml"
    shell: "bwa mem -5SP -T0 -t{threads} {input.fasta} {input.reads} -o {output}"
    #wrapper:
    #    "v1.1.0/bio/bwa/mem"


#rule pairtools_parse:
#    input: 
#        sam = "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.aligned.sam",
#        genome = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.genome"
#    output: 
#        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq40.aligned.pairsam"
#    log: 
#        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq40.aligned.log"
#    conda: "../envs/omni-c.yaml"
#    threads: 32
#    shell: "pairtools parse "
#           "--min-mapq 40 "
#           "--walks-policy 5unique "
#           "--max-inter-align-gap 30 "
#           "--nproc-in {threads} "
#           "--nproc-out {threads} "
#           "--chroms-path {input.genome} "
#           "{input.sam} > {output} 2> {log}"


rule pairtools_parse_mapq:
    input: 
        sam = "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.aligned.sam",
        genome = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.genome"
    output: 
        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.aligned.pairsam")
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.aligned.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: "pairtools parse "
           "--min-mapq {wildcards.mapq} "
           "--walks-policy 5unique "
           "--max-inter-align-gap 30 "
           "--nproc-in {threads} "
           "--nproc-out {threads} "
           "--chroms-path {input.genome} "
           "{input.sam} > {output} 2> {log}"


rule pairtools_sort:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.aligned.pairsam"
    output: 
        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.pairsam")
    log: 
        "logs/pairtools_sort/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        tmpdir=$(mktemp -d -p .)
        pairtools sort \
          --tmpdir=$tmpdir \
          --nproc {threads} \
          {input} > {output} 2> {log}
        rm -r $tmpdir
        """


rule pairtools_dedup:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.pairsam"
    output: 
        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.dedup.pairsam")
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.dedup.log"
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
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.dedup.pairsam"
    output: 
        bam = temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.bam"),
        pairs = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.pairs",
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.split.log"
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
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.bam"
    output:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam"
    params:
        #extra = "-m 100G"
    threads:  # Samtools takes additional threads through its option -@
        32     # This value - 1 will be sent to -@.
    wrapper:
        "v1.0.0/bio/samtools/sort"


rule samtools_index:
    input:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam"
    output:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam.bai"
    log:
        "output/Omni-C_samtools_index/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        32     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"
