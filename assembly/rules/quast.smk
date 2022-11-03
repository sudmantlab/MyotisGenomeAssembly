#configfile: "config.yaml"

def get_quast_inputs(wildcards, ccs_opts=config['ccs_opts']):
    hifi_path = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    input_samples = []
    for settings in ccs_opts:
        input_samples += [hifi_path.format(species=s[0], settings = settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples

rule quast:
    input:
        assemblies = expand("output/hifiasm-fasta/{species}/{settings}/{hifi_opts}/{species}.p_ctg.{hic}.fa", 
                            hifi_opts = config['hifi_opts'], settings = config['ccs_opts'], 
                            hic=config['conf_hic'], allow_missing=True)+ 
                     expand("output/yahs/{species}/{settings}/{hifi_opts}/{species}.p_ctg.{hic}_scaffolds_final.fa",
                            hifi_opts = config['hifi_opts'], settings = config['ccs_opts'],
                            hic=config['conf_hic'], allow_missing=True),
        pacbio = get_quast_inputs
    output: "output/quast/{species}/report.txt"
    conda: "../envs/genomeQC.yml"
    params:
        kmerSize = 21,
        pacbio = lambda wildcards: ["--pacbio {}".format(i) for i in get_quast_inputs(wildcards)],
        output_path = "output/quast/{species}/"
    shell: """ 
            quast.py \\
                {input.assemblies} \
                --eukaryote \\
                --large \\
                --k-mer-stats \\
                --k-mer-size {params.kmerSize} \\
                --rna-finding \\
                --conserved-genes-finding \\
                {params.pacbio} \\
                -t {threads} \\
                -o {params.output_path}
            """

