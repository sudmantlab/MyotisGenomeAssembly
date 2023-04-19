# bedtools sort -g data/chrom/mMyoLuc1.chrom data/gene_annotations/mMyoLuc1_RecSearch_1to1.bed > data/gene_annotations/mMyoLuc1_RecSearch_1to1.bed.sorted
# bedtools closest -D ref -a data/gene_annotations/mMyoLuc1_RecSearch_1to1.bed.sorted -b data/TEs/M_lucifugus_RM_full.bed.sorted > output/bedtools-closest/closest/mMyoLuc1-RMFull.bed
# bedtools shuffle -seed 452345 -i data/TEs/M_lucifugus_RM_full.gff.sorted -g data/chrom/mMyoLuc1.chrom | bedtools sort > data/TEs/M_lucifugus_RM_full.shuffled.sorted.bed
# bedtools closest -D ref -a data/gene_annotations/mMyoLuc1_RecSearch_1to1.bed.sorted -b data/TEs/M_lucifugus_RM_full.shuffled.sorted.bed> output/bedtools-closest/closest/mMyoLuc1-RMFull-shuffled.bed

rule bedtools_sort:
    input: 
        chrom = "data/chrom/{genome}.chrom",
        bed = "data/{class}/{genome}_{annotation}.bed"
    output: "data/{class}/{genome}_{annotation}.sorted.bed"
    conda: "../envs/bedtools.yaml"
    shell: "bedtools sort -g {input.chrom} -i {input.bed} > {output}"


rule bedtools_closest_default:
    input: 
        chrom = "data/chrom/{genome}.chrom",
        bedA = "data/gene_annotations/{genome}_{annotationA}.sorted.bed",
        bedB = "data/TEs/{genome}_{annotationB}.sorted.bed"
    conda: "../envs/bedtools.yaml"
    output: "output/bedtools-closest/default/{genome}-{annotationA}-{annotationB}.bed"
    shell: "bedtools closest -g {input.chrom} -D ref -io -a {input.bedA} -b {input.bedB} > {output}"


rule bedtools_shuffle:
    input: 
        chrom = "data/chrom/{genome}.chrom",
        bed = "data/{class}/{genome}_{annotation}.sorted.bed"
    output: "output/bedtools-shuffle/{class}/{seed}/{genome}_{annotation}.shuffled.bed"
    conda: "../envs/bedtools.yaml"
    shell: "bedtools shuffle -seed {wildcards.seed} -i {input.bed} -g {input.chrom} | bedtools sort -g {input.chrom} > {output}"

rule bedtools_closest_default_shuffleTEs:
    input: 
        chrom = "data/chrom/{genome}.chrom",
        bedA = "data/gene_annotations/{genome}_{annotationA}.sorted.bed",
        bedB = "output/bedtools-shuffle/TEs/{seed}/{genome}_{annotationB}.shuffled.bed"
    conda: "../envs/bedtools.yaml"
    output: "output/bedtools-closest/shuffled/default/{genome}-{annotationA}-shuffled{seed}_{annotationB}.bed"
    shell: "bedtools closest -D ref -io -g {input.chrom} -a {input.bedA} -b {input.bedB} > {output}"


rule bedtools_closest_default_shuffleGenes:
    input:
        chrom = "data/chrom/{genome}.chrom",
        bedA = "output/bedtools-shuffle/gene_annotations/{seed}/{genome}_{annotationA}.shuffled.bed",
        bedB = "data/TEs/{genome}_{annotationB}.sorted.bed"
    conda: "../envs/bedtools.yaml"
    output: "output/bedtools-closest/shuffled/default/{genome}-shuffled{seed}_{annotationA}-{annotationB}.bed"
    shell: "bedtools closest -D ref -io -g {input.chrom} -a {input.bedA} -b {input.bedB} > {output}"

rule bedtools_closest_default_shuffleAll:
    input:
        chrom = "data/chrom/{genome}.chrom",
        bedA = "output/bedtools-shuffle/gene_annotations/{seedA}/{genome}_{annotationA}.shuffled.bed",
        bedB = "output/bedtools-shuffle/TEs/{seedB}/{genome}_{annotationB}.shuffled.bed"
    conda: "../envs/bedtools.yaml"
    output: "output/bedtools-closest/shuffled/default/{genome}-shuffled{seedA}_{annotationA}-shuffled{seedB}_{annotationB}.bed"
    shell: "bedtools closest -D ref -io -g {input.chrom} -a {input.bedA} -b {input.bedB} > {output}"
