rule contigs_repeat_over_90pc:
    input:
        gff = '../data/repeatMasker/denovo/{genome}.denovo.out.reformated.gff',
        genome = 'data/genomes/{genome}.genome'
    output: 'output/contigs_to_drop/repeatMasker_denovo/{genome}.repeat_over_90pc'
    script: '../code/contigs_repeat_over_90pc.py'
