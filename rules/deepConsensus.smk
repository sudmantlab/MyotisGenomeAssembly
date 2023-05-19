rule ccs_all_minPasses_minRQ:
    input: "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{species}/all_minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs.bam",
        report="output/HiFi-CCS/{species}/all_minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs_report.txt"
    log: "logs/HiFi-deepConsensus/all_minPasses{minPasses}_minRQ{minRQ}/{species}/{pacbio1}/{pacbio2}/{id}.logs"
    threads: 52
    conda: "../envs/HiFiAssembly.yml"
    #singularity: "docker://google/deepconsensus:0.2.0-gpu"
    shell: "ccs --all -j {threads} --min-passes {wildcards.minPasses} --min-rq {wildcards.minRQ} --report-file {output.report} --log-file {log} {input} {output.ccs}"


rule subreads_to_ccs:
    input: 
        raw = "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam",
        ccs = "output/HiFi-CCS/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}.ccs.bam"
    output: 
        bam = "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}.subreads_to_ccs.bam",
        fasta = "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}.subreads_to_ccs.fasta"
    threads: 52
    conda: "../envs/HiFiAssembly.yml"
    #singularity: "docker://google/deepconsensus:0.2.0-gpu"
    shell: "actc -j {threads} "
           " {input.raw} "
           " {input.ccs} "
           " {output.bam}"

rule css_to_fasta:
    input: "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}.subreads_to_ccs.fasta"
    output: 
        fasta = "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}_ccs.fasta",
        fai = "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}_ccs.fasta.fai"
    threads: 52
    conda: "../envs/HiFiAssembly.yml"
    #singularity: "docker://google/deepconsensus:0.2.0-gpu"
    shell: """
        mv {input} {output.fasta}
        samtools faidx {output.fasta}
        """

rule deepConsensus:
    input:
        subreads_to_ccs = "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}.subreads_to_ccs.bam",
        ccs_fasta = "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}_ccs.fasta",
        model = "data/HiFi-deepConsensus/params.json"
    output:
        reads = "output/HiFi-deepConsensus/{species}/all_{ccs_opts}/{pacbio1}/{pacbio2}/{id}_ccs.fq"
    params:
        model = "data/HiFi-deepConsensus/checkpoint-50",
    threads: 52
    singularity: "docker://google/deepconsensus:0.2.0-gpu"
    shell: " deepconsensus run"
           " --subreads_to_ccs={input.subreads_to_ccs}"
           " --ccs_fasta={input.ccs_fasta}"
           " --checkpoint={params.model}"
           " --output={output.reads}"
           " --batch_zmws=100"


def get_hifiasm_inputs_DC(wildcards):
    hifi_path = "output/HiFi-deepConsensus/{species}/{settings}/{pacbio1}/{pacbio2}/{id}_ccs.fq"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    input_samples = [hifi_path.format(species=s[0], settings = wildcards.settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples


def get_hifiasm_hic_inputs_DC(wildcards):
    file_path = path_trimmed = "output/trimmed-hic/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    input_dict = {"left": [], "right": [], "hifi": get_hifiasm_inputs_DC(wildcards)}
    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
    samples = samples[samples.species==wildcards.species][samples.type == "Hi-C"]
    #print(samples)
    samples_grouped = samples.groupby(samples.sample_name)
    for sample in set(samples["sample_name"].tolist()):
        sample_subset = samples_grouped.get_group(sample)
        if len(sample_subset) == 0:
           raise Exception("No files available for sample {}".format(sample))
        input_dict["left"].extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R1"].itertuples()])
        input_dict["right"].extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R2"].itertuples()])
    return input_dict

rule hifiasm_l2_HiC_DC:
    version: "0.1"
    input: unpack(get_hifiasm_hic_inputs_DC)
    output:
        p_ctg_hap1 = multiext("output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.hic.hap1.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_ctg_hap2 = multiext("output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.hic.hap2.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_ctg = multiext("output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.hic.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_utg = multiext("output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.hic.p_utg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        r_utg = multiext("output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.hic.r_utg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        ec = "output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.ec.bin",
        # lk_bin = "output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.hic.lk.bin",
        ovlp_reverse = "output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm"
    threads: 52
    log: "output/hifiasm-HiC-DC/{species}/{settings}/l2/{species}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "hifiasm -o {params.prefix} --h1 {input.left} --h2 {input.right} -t {threads} -l2 {input.hifi} > {log} 2>&1"
