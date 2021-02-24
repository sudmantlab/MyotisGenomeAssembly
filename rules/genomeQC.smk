rule meryl_run:
    input: "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs.fastq.gz"
    output: directory("output/meryl/{species}/{pacbio1}/{pacbio2}/{id}.meryl")
    log: "logs/meryl/{species}/{pacbio1}/{pacbio2}/{id}.log"
    threads: 32
    shell: "meryl count threads={threads} k=21 {input} output {output} &> {log}"

rule merqury_run:
    input: 
        meryl = "output/meryl/{species}/{pacbio1}/{pacbio2}/{id}.meryl",
        p_ctg = "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm.p_ctg.gfa",
        a_ctg = "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm.a_ctg.gfa"        
    output: directory("output/merqury/{species}/{pacbio1}/{pacbio2}/{id}.asm")
    threads: 32
    shadow: "shallow"
    shell: "merqury.sh {input.meryl} {input.p_ctg} {input.a_ctg} {output}"

rule jellyfish_run:
    input: "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs.fastq.gz"
    output: "output/jellyfish/{species}/{pacbio1}/{pacbio2}/{id}.jf"
    threads: 32
    params:
        kmer = "21",
        memory = "1000000000"
    shell: "zcat {input} | jellyfish count -C -m {params.kmer} -s {params.memory} -t {threads} /dev/stdin -o {output}"

rule jellyfish_histo:
    input: "output/jellyfish/{species}/{pacbio1}/{pacbio2}/{id}.jf"
    output: "output/jellyfish/{species}/{pacbio1}/{pacbio2}/{id}.histo"
    threads: 32
    shell: "jellyfish histo -t {threads} {input} > {output}"
