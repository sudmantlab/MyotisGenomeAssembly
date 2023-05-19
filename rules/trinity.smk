def get_trinity_inputs(wildcards):
    path_trimmed = "output/trimmed/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    input_dict = {"left": [], "right": []}
    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
    samples = samples[samples.species==wildcards.species & samples.sample_name == wildcards.sample]
    print(samples)
    samples_grouped = samples.groupby(samples.sample_name)
    for sample in set(samples["sample_name"].tolist()):
        sample_subset = samples_grouped.get_group(sample)
        print(sample_subset)
        if len(sample_subset) == 0:
           raise Exception("No files available for sample {}".format(sample))
        input_dict["left"].extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R1"].itertuples()])
        input_dict["right"].extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R2"].itertuples()])
    return input_dict

rule trinity:
    input: 
        # unpack(get_trinity_inputs)
        left="output/trimmed/{species}/{sample}-trimmed_R1.fastq.gz",
        right="output/trimmed/{species}/{sample}-trimmed_R2.fastq.gz"
    output:
        "output/trinity/{species}/{sample}-trinity/Trinity.fasta"
    log:
        'logs/trinity/{species}/{sample}-trinity.log'
    params:
        extra="--SS_lib_type=RF"
    threads: 40
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_gb=192
    wrapper:
        "0.77.0/bio/trinity"


rule trinityStats: 
    input: "output/trinity/{species}/{sample}-trinity/Trinity.fasta"
    output: "output/trinity/{species}/{sample}-trinity/TrinityStats.txt"
    shell: "TrinityStats.pl {input} > {output}"

