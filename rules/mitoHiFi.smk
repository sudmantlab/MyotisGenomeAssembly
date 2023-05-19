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
    	genome = "data/genomes/{genome}.fa",
    	mt_fa = 'data/mitoHiFi/{genome}.fasta', #get_mitoHiFi_mitogenome_fasta,
    	mt_gb = 'data/mitoHiFi/{genome}.gb' #get_mitoHiFi_mitogenome_gb    	
    output: 
        stats =                "output/mitoHiFi/fromContig/{genome}/contigs_stats.tsv",
        final_genome_fasta =   "output/mitoHiFi/fromContig/{genome}/final_mitogenome.fasta",
        final_genome_genbank = "output/mitoHiFi/fromContig/{genome}/final_mitogenome.gb"
    params:
        output_dir = "output/mitoHiFi/fromContig/{genome}/",
        pct_identity = 90,
        circSize = 20000,
        circOffset = 1
    threads: 32
    shell: """
        export ROOTPROJDIR="$(pwd -P)"
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        singularity exec --bind $ROOTPROJDIR/:/lustre/ $ROOTPROJDIR/mitohifi_master.sif mitohifi.py \\
         -c  $ROOTPROJDIR/{input.genome} \\
         -f  $ROOTPROJDIR/{input.mt_fa} \\
         -g  $ROOTPROJDIR/{input.mt_gb} \\
         -t {threads} \\
         -a animal \\
         --circular-size {params.circSize} \\
         --circular-offset {params.circOffset} \\
         -p {params.pct_identity} \\
         -o 1
        cd $ROOTPROJDIR
    """
