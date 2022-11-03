def get_mitoHiFi_mitogenome_acc(wildcards):
    samples = pd.read_table("pep_mitochondrion.tsv", index_col="species")
    return samples.loc[wildcards.species].acc

def get_mitoHiFi_mitogenome_fasta(wildcards):
    return "data/mitoHiFi/" + get_mitoHiFi_mitogenome_acc(wildcards) + ".fasta"

def get_mitoHiFi_mitogenome_gb(wildcards):
    return "data/mitoHiFi/" + get_mitoHiFi_mitogenome_acc(wildcards) + ".gb"


rule mitoHiFi_reads:
    input: 
    	reads = get_hifiasm_inputs,
    	mt_fa = get_mitoHiFi_mitogenome_fasta,
    	mt_gb = get_mitoHiFi_mitogenome_gb    	
    output: 
        stats = "output/mitoHiFi/{species}/by_reads/{settings}/contigs_stats.tsv",
        final_genome_log = "output/mitoHiFi/{species}/by_reads/{settings}/final_mitogenome.annotation_MitoFinder.log",
        final_genome_fasta = "output/mitoHiFi/{species}/by_reads/{settings}/final_mitogenome.fasta",
        final_genome_genbank = "output/mitoHiFi/{species}/by_reads/{settings}/final_mitogenome.gb"
    #log:
    params:
        input = lambda wildcards: " ".join("$ROOTPROJDIR/{}".format(i) for i in get_hifiasm_inputs(wildcards)),
        output_dir = lambda wildcards: "$ROOTPROJDIR/output/mitoHiFi/{species}/by_reads/{settings}/".format(species=wildcards.species,settings=wildcards.settings),
        pct_identity = 70,
        circSize = 20000,
        circOffset = 3
    threads: 52
    singularity: "docker://docmanny/mitohifi:c06ed3e"
    shell: """
        ROOTPROJDIR="$(pwd -P)"
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        cat {params.input} > tmp_reads.fq.gz
        python $ROOTPROJDIR/code/MitoHiFi/mitohifi_v2.py \\
         -r tmp_reads.fq.gz \\
         -f  $ROOTPROJDIR/{input.mt_fa} \\
         -g  $ROOTPROJDIR/{input.mt_gb} \\
         -t {threads} \\
         -o 1 \\
         --circular-size {params.circSize} \\
         --circular-offset {params.circOffset} \\
         -p {params.pct_identity} \\
         -o 1
        rm tmp_reads.fq.gz
        cd $ROOTPROJDIR
    """


rule mitoHiFi_fromContig:
    input: 
    	genome = "output/hifiasm-fasta/{species}/{settings}/{opt}/{species}.p_ctg.hic.fa",
    	mt_fa = get_mitoHiFi_mitogenome_fasta,
    	mt_gb = get_mitoHiFi_mitogenome_gb    	
    output: 
        stats = "output/mitoHiFi/{species}/fromContig/{settings}/{opt}/contigs_stats.tsv",
        final_genome_log = "output/mitoHiFi/{species}/fromContig/{settings}/{opt}/final_mitogenome.annotation_MitoFinder.log",
        final_genome_fasta = "output/mitoHiFi/{species}/fromContig/{settings}/{opt}/final_mitogenome.fasta",
        final_genome_genbank = "output/mitoHiFi/{species}/fromContig/{settings}/{opt}/final_mitogenome.gb"
    #log: 
    params:
        output_dir = lambda wildcards: "output/mitoHiFi/{species}/fromContig/{settings}/{opt}/".format(species=wildcards.species,settings=wildcards.settings,opt=wildcards.opt),
        pct_identity = 80,
        circSize = 20000,
        circOffset = 1
    threads: 32
    singularity: "docker://docmanny/mitohifi:c06ed3e"
    shell: """
        ROOTPROJDIR="$(pwd -P)"
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        python $ROOTPROJDIR/code/MitoHiFi/mitohifi_v2.py \\
         -c  $ROOTPROJDIR/{input.genome} \\
         -f  $ROOTPROJDIR/{input.mt_fa} \\
         -g  $ROOTPROJDIR/{input.mt_gb} \\
         -t {threads} \\
         --circular-size {params.circSize} \\
         --circular-offset {params.circOffset} \\
         -p {params.pct_identity} \\
         -o 1
        cd $ROOTPROJDIR
    """
