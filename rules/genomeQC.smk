def get_meryl_union_sum_inputs(wildcards):
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    return ["output/meryl/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.meryl".format(species=s[0],
            settings = wildcards.settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]


def get_jellyfish_inputs(wildcards):
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    returnlist = ["output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz".format(species=s[0],
                  settings = wildcards.settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    #print(returnlist)
    return(returnlist)


rule meryl_run:
    input: "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz"
    output: directory("output/meryl/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.meryl")
    log: "logs/meryl/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.log"
    conda: "../envs/genomeQC.yml"
    threads: 52
    shell: "meryl count threads={threads} k=21 {input} output {output} &> {log}"

rule meryl_union_sum:
    input: get_meryl_union_sum_inputs
    output: directory("output/meryl/union-sum/{species}/{settings}/{species}.meryl")
    log: "logs/meryl-union-sum/{species}/{settings}/{species}.log"
    conda: "../envs/genomeQC.yml"
    threads: 52
    shell: "meryl union-sum threads={threads} {input} output {output} &> {log}"
    
rule merqury_run_haplotig:
    input: 
        meryl = "output/meryl/union-sum/{species}/{ccs_settings}/{species}.meryl",
        hap1 = "data/assemblies/{ccs_settings}/{hifiasm_opt}/{species}.hap1.{contigtype}.{hic}.fa",
        hap2 = "data/assemblies/{ccs_settings}/{hifiasm_opt}/{species}.hap2.{contigtype}.{hic}.fa"        
    output: directory("output/merqury/{ccs_settings}/{hifiasm_opt}/{species}.{contigtype}.{hic}.haplotigs")
    conda: "../envs/genomeQC.yml"
    threads: 52
    # shadow: "shallow"
    shell: """
        ROOTPROJDIR="$(pwd -P)"
        echo $ROOTPROJDIR
        mkdir -p {output}
        cd {output}
        pwd -P
        #rmdir $ROOTPROJDIR/{output}
        merqury.sh $ROOTPROJDIR/{input.meryl} $ROOTPROJDIR/{input.hap1} $ROOTPROJDIR/{input.hap2} {wildcards.species}.{wildcards.hic}.haplotig
        cd -
    """

rule merqury_run_yahs:
    input: 
        meryl = "output/meryl/union-sum/{species}/{ccs_settings}/{species}.meryl",
        consensus = "data/assemblies/{ccs_settings}/{hifiasm_opt}/{species}.{genometype}.{hic}_scaffolds_final.fa",
    output: directory("output/merqury/{species}/{ccs_settings}/{hifiasm_opt}/{species}.{hic}.{genometype}.yahs")
    conda: "../envs/genomeQC.yml"
    threads: 52
    # shadow: "shallow"
    shell: """
        ROOTPROJDIR="$(pwd -P)"
        echo $ROOTPROJDIR
        mkdir -p {output}
        cd {output}
        pwd -P
        #rmdir $ROOTPROJDIR/{output}
        merqury.sh $ROOTPROJDIR/{input.meryl} $ROOTPROJDIR/{input.consensus} {wildcards.species}.{wildcards.hic}.{wildcards.genometype}.yahs
        cd -
    """


rule merqury_run_arbitrary:
    input: 
        meryl = "output/meryl/union-sum/{species}/{ccs_settings}/{species}.meryl",
        consensus = "data/assemblies/{ccs_settings}/{hifiasm_opt}/{species}.{genometype}.{hic}.fa",
    output: directory("output/merqury/{ccs_settings}/{hifiasm_opt}/{species}.{genometype}.{hic}")
    conda: "../envs/genomeQC.yml"
    threads: 52
    # shadow: "shallow"
    shell: """
        ROOTPROJDIR="$(pwd -P)"
        echo $ROOTPROJDIR
        mkdir -p {output}
        cd {output}
        pwd -P
        #rmdir $ROOTPROJDIR/{output}
        merqury.sh $ROOTPROJDIR/{input.meryl} $ROOTPROJDIR/{input.consensus} {wildcards.species}.{wildcards.genometype}.{wildcards.hic}
        cd -
    """


rule jellyfish_run:
    input: get_jellyfish_inputs
    output: 
        bc = "output/jellyfish/{species}/{settings}/{species}.bc",
        jf = "output/jellyfish/{species}/{settings}/{species}.jf"
    conda: "../envs/genomeQC.yml"
    version: "0.0.2"
    threads: 52
    params:
        kmer = "21",
        generatorfile = lambda wildcards: "generators_{a}_{b}".format(a=wildcards.species, b=wildcards.settings)
    shell: """
        ls {input} | xargs -I[] echo zcat [] > {params.generatorfile}
        jellyfish bc -m {params.kmer} -s 100G -g {params.generatorfile} -G {threads} -t {threads} -o {output.bc} 
        jellyfish count -C -m {params.kmer} -s 3G -g {params.generatorfile} -G {threads} -t {threads} -o {output.jf} 
        rm {params.generatorfile}
    """

rule jellyfish_histo:
    input: "output/jellyfish/{species}/{settings}/{species}.jf"
    output: "output/jellyfish/{species}/{settings}/{species}.histo"
    conda: "../envs/genomeQC.yml"
    threads: 52
    shell: "jellyfish histo -t {threads} {input} > {output}"


rule genomescope:
    input: "output/jellyfish/{species}/{settings}/{species}.histo"
    output: 
        linearplot = "output/genomescope2/{species}/{settings}/{species}/linear_plot.png",
        logplot = "output/genomescope2/{species}/{settings}/{species}/log_plot.png",
        model = "output/genomescope2/{species}/{settings}/{species}/model.txt",
        progress = "output/genomescope2/{species}/{settings}/{species}/progress.txt",
        summary = "output/genomescope2/{species}/{settings}/{species}/summary.txt",
        translinearplot = "output/genomescope2/{species}/{settings}/{species}/transformed_linear_plot.png",
        translogplot = "output/genomescope2/{species}/{settings}/{species}/transformed_log_plot.png"
    params:
       prefix = "output/genomescope2/{species}/{settings}/{species}",
       kmer = 21
    log: "output/genomescope2/{species}/{settings}/{species}/log.txt"
    conda: "../envs/genomeQC.yml"
    shell: "genomescope2 -i {input} -o {params.prefix} -k {params.kmer} > {log} 2>&1"

def get_genomesize(gs2summary):
    import re
    genomesize = 0 
    pat = re.compile("(?<=len:)\d+")
    with open(gs2summary) as infile: 
        for line in infile:
            if pat.findall(line):
                genomesize = pat.search(line).group(0) 
    if genomesize: 
        return genomesize
    else:
        raise Exception("no genomesize found - did the file output from genomescope2 change?")

rule genometools:
    input: 
        fasta = "data/assemblies/{settings}/{opt}/{species}.{genometype}.{hic}.fa",
        gs2summary = "output/genomescope2/{species}/{settings}/{species}/log.txt"
    output: "output/gt-seqstat/{settings}/{opt}/{species}.{genometype}.{hic}.stats"
    run: 
        shell("./code/bin/gt "
              "seqstat -contigs -genome "
              "{genomesize} "
              "{input} > "
              "{output}".format(input=input.fasta,
                                output=output, 
                                genomesize=get_genomesize(input.gs2summary))
              )

rule genometools_yahs:
    input: 
        fasta = "output/yahs/{species}/{ccs_opts}/{hifiasm_opt}/{species}.{genometype}.{hic}_scaffolds_final.fa",
        gs2summary = "output/genomescope2/{species}/{ccs_opts}/{species}/log.txt"
    output: "output/gt-seqstat/{ccs_opts}/{hifiasm_opt}/{species}.{genometype}.{hic}.scaffold.stats"
    run: 
        shell("./code/bin/gt "
              "seqstat -contigs -genome "
              "{genomesize} "
              "-nstats 10 20 30 40 50 60 70 80 90 100 -- "
              "{input} > "
              "{output}".format(input=input.fasta,
                                output=output, 
                                genomesize=get_genomesize(input.gs2summary))
              )

#rule genometools_LJA:
#    input: 
#        fasta = "output/LJA/{species}/{ccs_opts}/{hifiasm_opt}/{species}.{genometype}.{hic}_scaffolds_final.fa",
#        gs2summary = "output/genomescope2/{species}/{ccs_opts}/{species}/log.txt"
#    output: "output/gt-seqstat/{species}/{ccs_opts}/{hifiasm_opt}/{species}.{genometype}.{hic}.scaffold.stats"
#    run: 
#        shell("./code/bin/gt "
#              "seqstat -contigs -genome "
#              "{genomesize} "
#              "{input} > "
#              "{output}".format(input=input.fasta,
#                                output=output, 
#                                genomesize=get_genomesize(input.gs2summary))
#              )

rule symlink_assemblies:
    input: "output/{assembler}/{species}/{ccs_settings}/{assembler_opts}/{species}.{assembly}.fa"
    output: "data/assemblies/{species}/{ccs_settings}/{assembler}-{assembler_opts}/{species}.{assembly}.fa"
    shell: "ln -sf '$(pwd -P)'/{input} {output}"


rule genometools_generalized:
    input: 
        fasta = "data/assemblies/{ccs_settings}/{hifiasm_opts}/{species}.{assembly}.fa",
        gs2summary = "output/genomescope2/{species}/{ccs_settings}/{species}/log.txt"
    output: "output/gt-seqstat/{ccs_settings}/{hifiasm_opts}/{species}.{assembly}.stats"
    run: 
        shell("./code/bin/gt "
              "seqstat -contigs -genome "
              "{genomesize} "
              "-nstats 10 50 90 -- "
              "{input} > "
              "{output}".format(input=input.fasta,
                                output=output, 
                                genomesize=get_genomesize(input.gs2summary))
              )
