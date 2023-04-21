rule iqtree_mfp:
    input:
        'data/genes/{class}/{gene}.fa'
    output:
        'output/iqtree/{class}/{gene}/{gene}.iqtree'
    params:
        prefix = 'output/iqtree/{class}/{gene}/{gene}',
        seqtype = 'CODON',
        opts = '-alrt 1000 --wbt -B 1000 -bnni'
    conda: "../envs/phylo.yaml"
    threads: 10
    shell: 
        'iqtree -safe -cptime 2 -s {input} -T {threads} -st {params.seqtype} -pre {params.prefix} {params.opts}'

#rule iqtree_
