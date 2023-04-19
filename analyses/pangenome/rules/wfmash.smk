rule wfmash_approx_1to1:
    version: 0.1
    input: 
        ref_genome = "data/genomes/{ref_genome}.fa",
        query_genome = "data/genomes/{query_genome}.fa"
    output: "output/wfmash/1-to-1/{ref_genome}_{query_genome}.approx.paf"
    threads: 32
    #conda: "../envs/wfmash.smk"
    #shell: "wfmash {input.ref_genome} {input.query_genome} -t40 -k21 -m > {output}"
    shell: "singularity run wfmash-docker_0.10.0.sif wfmash {input.ref_genome} {input.query_genome} -t {threads} -k21 -m > {output}"

rule wfmash_1to1_refine:
    version: 0.1
    input: 
        paf = "output/wfmash/1-to-1/{ref_genome}_{query_genome}.approx.paf",
        ref_genome = "data/genomes/{ref_genome}.fa",
        query_genome = "data/genomes/{query_genome}.fa"
    output: "output/wfmash/1-to-1/{ref_genome}_{query_genome}.paf"
    threads: 30
    #conda: "../envs/wfmash.smk"
    #shell: "wfmash {input.ref_genome} {input.query_genome} -t40 -k21 -i {input.paf} > {output}"
    shell: "singularity run wfmash-docker_0.10.0.sif wfmash {input.ref_genome} {input.query_genome} -t40 -k21 -i {input.paf} > {output}"
    

genomes = ["mMyoAui1", "mMyoCai1", "mMyoEvo1", "mMyoLuc1", "mMyoOcc1", "mMyoThy1", "mMyoVel1", "mMyoVol1", "mMyoYum1", "mMyoSep1"]
genomes_noYum = ["mMyoAui1", "mMyoCai1", "mMyoEvo1", "mMyoLuc1", "mMyoOcc1", "mMyoThy1", "mMyoVel1", "mMyoVol1"]
haps = ["hap1","hap2"]

rule wfmash_approx_AvA_refOnly:
    version: 0.1
    input: 
        expand("data/genomes/{genome}.cleaned.hapheader.fa", genome=genomes)
    output: 
        "output/wfmash/AvA/neartic_myotis_refOnly.approx.paf"
    threads: 32
    params:
        n_map = len(genomes) - 1
    #singularity: "docker://docmanny/wfmash-docker:0.10.0"
    #shell: "wfmash {input} -t {threads} -X -n {params.n_map} -m > {output}"
    shell: "singularity run wfmash-docker_0.10.0.sif wfmash {input} -t {threads} -X -n {params.n_map} -m > {output}"


rule wfmash_approx_AvA_refAndHap:
    version: 0.1
    input: 
        expand("data/genomes/{genome}.cleaned.hapheader.fa", genome=genomes) + 
        expand("data/genomes/{genome}.{hap}.chr_M.hapheader.fa", genome=genomes_noYum, hap=haps) + 
        expand("data/genomes/mMyoYum1.{hap}.hapheader.fa", hap=haps)
    output: 
        "output/wfmash/AvA/neartic_myotis_refAndHap.approx.paf"
    threads: 32
    #singularity: "docker://docmanny/wfmash-docker:0.10.0"
    #shell: "/opt/conda/bin/wfmash {input} -t{threads} -k21 -m > {output}"
    shell: "singularity run wfmash-docker_0.10.0.sif wfmash {input} -t{threads} -k21 -m > {output}"


rule mash_triangle_scaffold:
    version: 0.1
    input:
        expand("data/genomes/{genome}.cleaned.hapheader/{genome}.cleaned.hapheader.{{contig}}.fa", genome=genomes_noYum)
    output: "output/mash-triangle/neartic_myotis.mash_triangle.{contig}.txt"
    threads: 32
    params:
        kmer = 19
    #conda: 
    shell: "mash triangle -p {threads} -k {params.kmer} {input} > {output}"

rule mash_triangle_genome:
    version: 0.1
    input:
        expand("data/genomes/{genome}.cleaned.hapheader.fa", genome=genomes)
    output: "output/mash-triangle/neartic_myotis.mash_triangle.all.txt"
    threads: 32
    params:
        kmer = 19
    #conda: 
    shell: "mash triangle -p {threads} -k {params.kmer} {input} > {output}"

