rule hardmask_repeatmasker_union:
    input: 
        fasta = "data/genomes/{genome}.cleaned.fa",
        rm = "data/RepeatMasker/union/{genome}.full.denovo.gff"
    output: "output/genomes_hardmask/{genome}.cleaned.hardmask.fa"
    shell: "bedtools maskfasta -fullHeader -fi {input.fasta} -fo {output} -bed {input.rm}"


rule sedef:
    input: "output/genomes_hardmask/{genome}.cleaned.hardmask.fa"
    output: 
        dir = directory("output/sedef/{genome}"),
        bed = "output/sedef/{genome}/final.bed"
    threads: 40
    shell: "./code/sedef/sedef.sh -o {output.dir} -j {threads} -f {input}"


rule get_mainChr_hardmask:
    input: "output/genomes_hardmask/{genome}.cleaned.hardmask.fa"
    output: "output/genomes_hardmask/{genome}.cleaned.hardmask.mainChr.fa"
    params:
        mainChr = " ".join(["SUPER__" + str(i) for i in range(1,24)])
    shell: "samtools faidx {input} {params.mainChr} > {output}"


