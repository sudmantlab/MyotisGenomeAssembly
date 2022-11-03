import os
import pandas as pd

def get_lja_inputs(wildcards):
    hifi_path = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    input_samples = [hifi_path.format(species=s[0], settings = wildcards.settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples

rule laJollaAssembler:
    version: "0.2"
    input: get_lja_inputs
    output: 
        fasta = "output/LJA/{species}/{ccs_settings}/k{k_err}_K{k_multiDBG}_diploid/assembly.fasta",
        gfa = "output/LJA/{species}/{ccs_settings}/k{k_err}_K{k_multiDBG}_diploid/mdbg.gfa",
        fasta_homo = "output/LJA/{species}/{ccs_settings}/k{k_err}_K{k_multiDBG}_diploid/mdbg/assembly.hpc.fasta",
        gfa_homo = "output/LJA/{species}/{ccs_settings}/k{k_err}_K{k_multiDBG}_diploid/mdbg/mdbg.hpc.gfa"
    params:
        prefix = "output/LJA/{species}/{ccs_settings}/k{k_err}_K{k_multiDBG}_diploid/",
        reads = lambda wildcards: ["--reads " + i for i in get_lja_inputs(wildcards)] 
    threads: 52
    log: "output/LJA/{species}/{settings}/dbg.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "code/LJA/bin/lja --diploid -o {params.prefix} -t {threads} -k {wildcards.k_err} -K {wildcards.k_multiDBG} {params.reads}"
