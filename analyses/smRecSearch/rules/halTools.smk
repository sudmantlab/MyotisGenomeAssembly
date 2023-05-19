# Functions
import pandas as pd


def genome_to_species_hal(genome, df):
    db = pd.read_csv(df)
    result = db[db.UCSCName == genome].species.to_string(index=False).strip()
    return result

# Rules

rule getScaffoldsWithDups:
    input: "output/finalBed/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_duplicatesOnly_GenBankAccn.bed"
    output: "output/chrWithDups/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_chr.txt"
    shell: "cut -f1 {input} | sort | uniq > {output}"


rule cleanGenomeTable:
    input: 
      hal_db = "data/hal/200m_species.tsv",
      ncbi_db = "data/hal/ncbi_db.csv"
    output: "data/hal/200m_species_clean.tsv"
    conda: "../envs/conda_cactus.yaml"
    script: "../code/cleanGenomeTable.py"
        

rule getSyntenyChr:
    input: 
        hal = "data/hal/200m-v1.hal",
        species_genome = "data/hal/200m_species_clean.tsv"
    output: temp("output/syntenyBlocks/psl/{genomeA}-{genomeB}-{chr}_maxAnchorDistance{maxAnchorDistance}_minBlockSize{minBlockSize}.psl")
    params:
        SpeciesA = lambda wildcards, input: genome_to_species_hal(wildcards.genomeA, input.species_genome),
        SpeciesB = lambda wildcards, input: genome_to_species_hal(wildcards.genomeB, input.species_genome)
    conda: "../envs/conda_cactus.yaml"
    shell: "halSynteny --queryGenome {params.SpeciesA} --targetGenome {params.SpeciesB} --maxAnchorDistance {wildcards.maxAnchorDistance} --minBlockSize {wildcards.minBlockSize} --queryChromosome {wildcards.chr} {input} {output}"


rule expandLocus:
    input: 
      bed= "output/finalBed/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}.bed",
      genome = "data/genome/{genome}.genome"
    output: "output/BED_expandedLocus/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_left{lshift}_right{rshift}.bed"
    conda: "../envs/conda_cactus.yaml"
    threads: 28
    shell: "cat {input.bed} | parallel -j {threads} 'echo {{}} | bedtools slop -g {input.genome} -l {wildcards.lshift} -r {wildcards.rshift} -s' > {output}"
    
    
rule mergeLocus:
    input: 
        bed = "output/BED_expandedLocus/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{hit_type}_{type}_left{lshift}_right{rshift}.bed",
        genome = "data/genome/{genome}.genome"
    output: "output/mergedLocus/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{hit_type}_{type}_left{lshift}_right{rshift}_mergedLocus.bed"
    conda: "../envs/conda_cactus.yaml"
    shell: "cat {input.bed} | bedtools sort -g {input.genome} | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > {output}"

rule UCSCScaffoldsToGenBank:
    input:
        bed = "output/mergedLocus/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}.bed",
        replace_table = "data/assemblyReports/{genome}_assemblyReport.tsv"
    output: "output/BED_GenBankAccn/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_GenBankAccn.bed"
    params:
        frm = "UCSC_style_name"
    conda: "../envs/conda_cactus.yaml"
    script: "../code/convert_ucsc_scaffold_to_GenBank.py"


# rule GenBanktoUCSCScaffolds:
#     input:
#         bed = "output/BED_GenBankAccn/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_GenBankAccn.bed",
#         replace_table = "data/assemblyReports/{genome}_assemblyReport.tsv"
#     output: "output/BED_UCSCAccn/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_USCSCAccn.bed"
#     params:
#           frm = "GenBankAccn"
#     conda: "../envs/conda_cactus.yaml"
#     script: "../code/convert_ucsc_scaffold_to_GenBank.py"


rule halLiftOver:
    input: 
        hal = "data/hal/200m-v1.hal",
        species_genome = "data/hal/200m_species_clean.tsv",
        bed = "output/BED_GenBankAccn/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}_{hit_type}_{type}_GenBankAccn.bed"
    output: "output/halLiftOver/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_{type}_GenBankAccn.bed"
    params:
        SpeciesA = lambda wildcards, input: genome_to_species_hal(wildcards.genomeA, input.species_genome),
        SpeciesB = lambda wildcards, input: genome_to_species_hal(wildcards.genomeB, input.species_genome)
    conda: "../envs/conda_cactus.yaml"
    threads: 28
    shell: "cat {input.bed} | parallel -j {threads} 'echo {{}} | halLiftover {input.hal} {params.SpeciesA} /dev/stdin {params.SpeciesB} /dev/stdout' > {output}"
    #shell: "halLiftover {input.hal} {params.SpeciesA} {input.bed} {params.SpeciesB} {output}"


rule mergeLiftOver:
    input: 
        bed = "output/halLiftOver/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_{type}_GenBankAccn.bed",
        genome = "data/genome/{genomeB}.genome"
    output: "output/halLiftOver/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_{type}_GenBankAccn_merged.bed"
    conda: "../envs/conda_cactus.yaml"
    shell: "cat {input.bed} | bedtools sort -g {input.genome} | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > {output}"


rule reciprocalLiftOver:
    input: 
        hal = "data/hal/200m-v1.hal",
        species_genome = "data/hal/200m_species_clean.tsv",
        bed = "output/halLiftOver/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_{type}_GenBankAccn_merged.bed"
    output: "output/recHalLiftOver/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}-to-{genomeA}_{hit_type}_{type}_GenBankAccn_merged.bed"
    params:
        SpeciesA = lambda wildcards, input: genome_to_species_hal(wildcards.genomeA, input.species_genome),
        SpeciesB = lambda wildcards, input: genome_to_species_hal(wildcards.genomeB, input.species_genome)
    conda: "../envs/conda_cactus.yaml"
    threads: 28
    shell: "cat {input.bed} | parallel -j {threads} 'echo {{}} | halLiftover {input.hal} {params.SpeciesB} /dev/stdin {params.SpeciesA} {output}"


rule fastaFromBed_mergedLocus:
    input:
        bed = "output/mergedLocus/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{hit_type}_{type}_left{lshift}_right{rshift}_mergedLocus.bed",
        fasta = "data/genome/{genome}.fa",
        fai = "data/genome/{genome}.fa.fai",
        genome = "data/genome/{genome}.genome"
    output: "output/duplicateGeneFASTA/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{hit_type}_{type}_left{lshift}_right{rshift}_mergedLocus.fa"
    conda: "../envs/conda_cactus.yaml"
    shell: "cat {input.bed} | bedtools sort -g {input.genome} | fastaFromBed -name+ -fullHeader -fi {input.fasta} -bed - > {output}"
    
rule fastaFromBed_exons:
    input:
        bed = "output/finalBed/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{hit_type}_duplicatesOnly.bed",
        fasta = "data/genome/{genome}.fa",
        fai = "data/genome/{genome}.fa.fai",
        genome = "data/genome/{genome}.genome"
    output: "output/duplicateGeneFASTA/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{hit_type}_duplicatesOnly.fa"
    conda: "../envs/conda_cactus.yaml"
    shell: "cat {input.bed} | bedtools sort -g {input.genome} | fastaFromBed -name+ -fullHeader -split -fi {input.fasta} -bed - > {output}"
