import os
import pandas as pd

rule ccs_consensus_noopts:
    input: "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{species}/defaults/{pacbio1}/{pacbio2}/{id}.ccs.bam",
        report="output/HiFi-CCS/{species}/defaults/{pacbio1}/{pacbio2}/{id}.ccs_report.txt"
    log: "logs/HiFi-CCS/{species}/defaults/{pacbio1}/{pacbio2}/{id}.logs"
    threads: 32
    conda: "../envs/HiFiAssembly.yml"
    shell: "ccs -j {threads} --report-file {output.report} --log-file {log} {input} {output.ccs}"


rule ccs_consensus_minPasses_minRQ:
    input: "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{species}/minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs.bam",
        report="output/HiFi-CCS/{species}/minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs_report.txt"
    log: "logs/HiFi-CCS/minPasses{minPasses}_minRQ{minRQ}/{species}/{pacbio1}/{pacbio2}/{id}.logs"
    threads: 32
    conda: "../envs/HiFiAssembly.yml"
    shell: "ccs -j {threads} --min-passes {wildcards.minPasses} --min-rq {wildcards.minRQ} --report-file {output.report} --log-file {log} {input} {output.ccs}"


rule HiFiAdapterFilt:
    input: "output/HiFi-CCS/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.bam"
    output: 
        #blocklist = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.blocklist", 
        #contaminant = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.contaminant.blastout", 
        #fasta = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.fasta", 
        #fastq = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.fastq",
        #filtered = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq"
        filtered = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz",
        stats = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.stats",
    log: "logs/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.log"
    version: 2.00
    params:
        outDir = lambda wildcards, output: os.path.dirname(output[0]),
        inBaseName = lambda wildcards, input: os.path.basename(input[0]),
        inDir = lambda wildcards, input: os.path.dirname(input[0]),
        inPref = lambda wildcards, input: os.path.splitext(os.path.basename(input[0]))[0],        
    conda: "../envs/HiFiAssembly.yml"
    threads: 10
    shell: """
        ROOTPROJDIR=$(pwd -P)
        cd {params.inDir}
        #mkdir $ROOTPROJDIR/{params.outDir}
        pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
        cd -
    """

rule HiFiAdapterFilt_fastq:
    input: "output/HiFi-CCS/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.fastq"
    output: 
        #blocklist = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.blocklist", 
        #contaminant = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.contaminant.blastout", 
        #fasta = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.fasta", 
        #fastq = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.fastq",
        #filtered = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq"
        filtered = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz",
        stats = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.stats",
    log: "logs/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.log"
    #version: 2.00
    params:
        outDir = lambda wildcards, output: os.path.dirname(output[0]),
        inBaseName = lambda wildcards, input: os.path.basename(input[0]),
        inDir = lambda wildcards, input: os.path.dirname(input[0]),
        inPref = lambda wildcards, input: os.path.splitext(os.path.basename(input[0]))[0],        
    conda: "../envs/HiFiAssembly.yml"
    threads: 10
    shell: """
        ROOTPROJDIR=$(pwd -P)
        cd {params.inDir}
        #mkdir $ROOTPROJDIR/{params.outDir}
        pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
        cd -
    """


rule bam2fastq:
    input: "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs.bam"
    output: "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs.fastq.gz"
    params:
        prefix = "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs"
    conda: "../envs/HiFiAssembly.yml"
    shell: "bam2fastq -o {params.prefix} {input}"


def get_hifiasm_inputs(wildcards):
    hifi_path = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    input_samples = [hifi_path.format(species=s[0], settings = wildcards.settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples

rule hifiasm_noopts:
    version: "0.2.2"
    input: get_hifiasm_inputs
    output: 
        p_ctg_hap1 = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.p_ctg.gfa",
        p_ctg_lowQ = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.p_ctg.noseq.gfa",
        p_utg = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.p_utg.gfa",
        p_utg_lowQ = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.p_utg.lowQ.bed",
        p_utg_noseq = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.p_utg.noseq.gfa",
        r_utg = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.r_utg.gfa",
        r_utg_lowQ = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.r_utg.lowQ.bed",
        r_utg_noseq = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.bp.r_utg.noseq.gfa",
        ec = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.ec.bin",
        ovlp_reverse = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm/{species}/{settings}/no_opts/{species}.asm"
    threads: 32
    log: "output/hifiasm/{species}/{settings}/no_opts/{species}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    # shell: "hifiasm -l2 -o {params.prefix} -t {threads} {input} > {log} 2>&1"
    shell: "hifiasm -o {params.prefix} -t {threads} {input} > {log} 2>&1"

rule hifiasm_l2:
    version: "0.2.2"
    input: get_hifiasm_inputs
    output: 
        p_ctg_hap1 = multiext("output/hifiasm/{species}/{settings}/l2/{species}.asm.bp.hap1.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_ctg_hap2 = multiext("output/hifiasm/{species}/{settings}/l2/{species}.asm.bp.hap2.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_ctg = multiext("output/hifiasm/{species}/{settings}/l2/{species}.asm.bp.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_utg = multiext("output/hifiasm/{species}/{settings}/l2/{species}.asm.bp.p_utg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        r_utg = multiext("output/hifiasm/{species}/{settings}/l2/{species}.asm.bp.r_utg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        ec = "output/hifiasm/{species}/{settings}/l2/{species}.asm.ec.bin",
        #lk_bin = "output/hifiasm/{species}/{settings}/l2/{species}.asm.hic.lk.bin",
        ovlp_reverse = "output/hifiasm/{species}/{settings}/l2/{species}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm/{species}/{settings}/l2/{species}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm/{species}/{settings}/l2/{species}.asm"
    threads: 32
    log: "output/hifiasm/{species}/{settings}/l2/{species}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "hifiasm -l2 -o {params.prefix} -t {threads} {input} > {log} 2>&1"


def get_hifiasm_hic_inputs(wildcards):
    file_path = path_trimmed = "output/trimmed-hic/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    input_dict = {"left": [], "right": [], "hifi": get_hifiasm_inputs(wildcards)}
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


rule hifiasm_noopts_HiC:
    version: "0.2.1"
    input: unpack(get_hifiasm_hic_inputs)
    output:
        p_ctg_hap1 = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.p_ctg.gfa",
        p_ctg_lowQ = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.p_ctg.noseq.gfa",
        p_utg = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.p_utg.gfa",
        p_utg_lowQ = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.p_utg.lowQ.bed",
        p_utg_noseq = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.p_utg.noseq.gfa",
        r_utg = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.r_utg.gfa",
        r_utg_lowQ = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.r_utg.lowQ.bed",
        r_utg_noseq = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.r_utg.noseq.gfa",
        ec = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.ec.bin",
        #lk_bin = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.hic.lk.bin",
        ovlp_reverse = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm"
    threads: 32
    log: "output/hifiasm-HiC/{species}/{settings}/no_opts/{species}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    # shell: "hifiasm -o {params.prefix} --h1 {input.left} --h2 {input.right} -t {threads} -l2 {input.hifi} > {log} 2>&1"
    shell: "hifiasm -o {params.prefix} --h1 {input.left} --h2 {input.right} -t {threads} {input.hifi} > {log} 2>&1"


rule hifiasm_l2_HiC:
    version: "0.2.2"
    input: unpack(get_hifiasm_hic_inputs)
    output:
        p_ctg_hap1 = multiext("output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.hic.hap1.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_ctg_hap2 = multiext("output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.hic.hap2.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_ctg = multiext("output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.hic.p_ctg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        p_utg = multiext("output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.hic.p_utg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        r_utg = multiext("output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.hic.r_utg", ".gfa", ".lowQ.bed", ".noseq.gfa"),
        ec = "output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.ec.bin",
        # lk_bin = "output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.hic.lk.bin",
        ovlp_reverse = "output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm"
    threads: 32
    log: "output/hifiasm-HiC/{species}/{settings}/l2/{species}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "hifiasm -o {params.prefix} --h1 {input.left} --h2 {input.right} -t {threads} -l2 {input.hifi} > {log} 2>&1"


rule gfaToFa:
    input: "output/hifiasm/{species}/{settings}/{opt}/{species}.asm.bp.{genometype}.gfa"
    output: 
        fa = "output/hifiasm-fasta/{species}/{settings}/{opt}/{species}.{genometype}.bp.fa",
        ln =  "data/assemblies/{settings}/{opt}/{species}.{genometype}.bp.fa"
    log: "logs/gfaToFa/{species}/{settings}/{opt}/{species}.{genometype}.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "gfatools gfa2fa {input} > {output.fa} 2> {log}; "
           "mkdir -p data/assemblies/{wildcards.settings}/{wildcards.opt}/; "
           "ln -s $(pwd -P)/{output.fa} {output.ln}"

rule gfaToFa_hic:
    input: "output/hifiasm-HiC/{species}/{settings}/{opt}/{species}.asm.hic.{genometype}.gfa"
    output: 
        fa = "output/hifiasm-fasta/{species}/{settings}/{opt}/{species}.{genometype}.hic.fa",
        ln =  "data/assemblies/{settings}/{opt}/{species}.{genometype}.hic.fa"
    log: "logs/gfaToFa-HiC/{species}/{settings}/{opt}/{species}.{genometype}.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "gfatools gfa2fa {input} > {output.fa} 2> {log}; "
           "mkdir -p data/assemblies/{wildcards.settings}/{wildcards.opt}/; "
           "ln -s $(pwd -P)/{output.fa} {output.ln}"


