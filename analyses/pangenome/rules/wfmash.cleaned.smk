rule wfmash_approx_1to1:
    version: 0.1
    input: 
        ref_genome = "data/genomes/{ref_genome}.cleaned.fa",
        query_genome = "data/genomes/{query_genome}.cleaned.fa"
    output: "output/wfmash/1-to-1/{ref_genome}.cleaned_{query_genome}.cleaned.approx.paf"
    threads: 32
    shell: "singularity run wfmash-docker_0.10.0.sif wfmash {input.ref_genome} {input.query_genome} -t {threads} -k21 -m > {output}"

rule wfmash_1to1_refine:
    version: 0.1
    input: 
        paf = "output/wfmash/1-to-1/{ref_genome}.cleaned_{query_genome}.cleaned.approx.paf",
        ref_genome = "data/genomes/{ref_genome}.cleaned.fa",
        query_genome = "data/genomes/{query_genome}.cleaned.fa"
    output: "output/wfmash/1-to-1/{ref_genome}.cleaned_{query_genome}.cleaned.paf"
    threads: 30
    shell: "singularity run wfmash-docker_0.10.0.sif wfmash {input.ref_genome} {input.query_genome} -t40 -k21 -i {input.paf} > {output}"

