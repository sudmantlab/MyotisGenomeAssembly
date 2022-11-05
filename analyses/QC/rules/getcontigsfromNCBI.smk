rule pipe_zcat:
    input: "data/genomes/{genome}.fa.gz"
    output: "data/genomes/{genome}.fa"
    shell: "zcat {input} > {output}"


rule samtools_faidx:
    input: "data/genomes/{genome}.fa"
    output: "data/genomes/{genome}.fa.fai"
    shell: "samtools faidx {input}"

rule genome_file:
    input: "data/genomes/{genome}.fa.fai"
    output: "data/genomes/{genome}.genome"
    shell: "cut -f1,2 {input} > {output}"

rule nbed:
    input: "data/2bit/{genome}.2bit"
    output: "data/nBed/{genome}_N.bed"
    shell: "code/twoBitInfo -nBed {input} {output}"

rule bedtools_complement:
    input: 
        chrom = "data/genomes/{genome}.genome",
        gap = "data/nBed/{genome}_N.bed"
    output: "output/genomes_contigs/{genome}.contigs.bed"
    shell: "bedtools complement -g {input.chrom} -i {input.gap} > {output}"

rule get_fasta_contigs:
    input: 
        fa = "data/genomes/{genome}.fa",
        fai = "data/genomes/{genome}.fa.fai",
        bed = "output/genomes_contigs/{genome}.contigs.bed"
    output: "output/genomes_contigs/{genome}.contigs.fa"
    shell: "bedtools getfasta -name+ -fullHeader -fi {input.fa} -bed {input.bed} -fo {output}"

rule link_fasta_contigs:
    input: "output/genomes_contigs/{genome}.contigs.fa"
    output: "data/genomes/{genome}.contigs.fa"
    shell: "PROJDIR=$(pwd -P); ln -sf $PROJDIR/{input} {output}"
