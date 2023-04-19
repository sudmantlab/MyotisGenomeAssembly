#ruleorder: gz_compress > gz_extract > faToTwoBit > twoBitToFa

#rule twoBitToFa:
#    input: ancient("data/2bit/{genome}.2bit")
#    output: protected("data/genome/{genome}.fa")
#    conda: "../envs/conda_twoBitToFa.yaml"
#    shell:
#        "twoBitToFa {input} {output}"

rule faToTwoBit:
    input: "data/genome/{genome}.fa"
    output: "data/2bit/{genome}.2bit"
    conda: "../envs/conda_faToTwoBit.yaml"
    shell:
        "faToTwoBit {input} {output}"

rule fasta_index:
    input: "data/genome/{genome}.fa"
    output: "data/genome/{genome}.fa.fai"
    conda: "../envs/conda_tuxedo.yaml"
    shell: "samtools faidx {input}"
    
rule genome_file:
    input: "data/genome/{genome}.fa.fai"
    output: "data/genome/{genome}.genome"
    shell: "awk -v OFS='\t' {{'print $1,$2'}} {input} > {output}"

#rule gz_compress:
#    input: "{file}.fa"
#    output: "{file}.fa.gz"
#    wildcard_constraints:
#        extension="[A-Za-z0-9]+(?!gz)"
#    shell: "gzip -k {input}"

rule manualBlat:
    input:
        twoBit="data/2bit/{genomeB}.2bit",
        query="output/duplicateGeneFASTA/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}_{hit_type}_duplicatesOnly.fa"
    output: "output/blat/forwardPSL/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_duplicatesOnly_filtered.psl"
    conda: "../envs/conda_blat.yaml"
    shell: "blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 {input.twoBit} {input.query} {output}"
    
rule PslToBED:
    input: "output/blat/forwardPSL/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_duplicatesOnly_filtered.psl"
    output: "output/blat/forwardBED/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_duplicatesOnly_filtered.bed"
    conda: "../envs/conda_blat.yaml"
    shell: "psl2bed-megarow < {input} > {output}"
    
rule BedToFasta:
    input: 
        bed = "output/blat/forwardBED/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_duplicatesOnly_filtered.bed",
        fasta = "data/genome/{genomeB}.fa",
        fai = "data/genome/{genomeB}.fa.fai",
        genome = "data/genome/{genomeB}.genome"
    output: "output/blat/forwardFasta/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}-to-{genomeB}_{hit_type}_duplicatesOnly_filtered.fa"
    conda: "../envs/conda_blat.yaml"
    shell: "cat {input.bed} | bedtools sort -g {input.genome} | fastaFromBed -name+ -fullHeader -split -fi {input.fasta} -bed - > {output}"


def getGeneCopyTableRBB(wildcards):
    try:
        search_path = "output/geneCopyTable/{searchType}".format(searchType=wildcards.searchType) + "/{genome}_geneCopyTable_RBB_filtered.tsv"
        RBBGenomes, = glob_wildcards(search_path)
        files = expand(search_path, genome = RBBGenomes)
        if files:
            return files
        else:
            raise Exception("No files found! Files: {0}\nSearch Path: {1}".format(str(files), search_path))
    except Exception:
        raise Exception("No files found! Files: {0}\nSearch Path: {1}".format(str(files), search_path))

def getGeneCopyTableRBB_Atlantogenata(wildcards):
    search_path = "output/geneCopyTable/{searchType}".format(searchType=wildcards.searchType) + "/{genome}_geneCopyTable_RBB_filtered.tsv"
    RBBGenomes = ["eleMaxD","loxAfr4","loxCycF","mamAmeI","mamColU","mamPriV","palAntN","proCap-Pcap_2.0_HiC","triManLat2","choHof-C_hoffmanni-2.0.1_HiC","chrAsi1m","dasNov3","eleEdw1m", "oryAfe2", "echTel2"]
    return expand(search_path, genome = RBBGenomes)
    
    
def getGeneCopyTableEvidenced_Atlantogenata(wildcards):
    search_path = "output/geneCopyTable/{searchType}".format(searchType=wildcards.searchType) + "/{genome}_geneCopyTable_evidenced_filtered.tsv"
    RBBGenomes = ["eleMaxD","loxAfr4","loxCycF","mamAmeI","mamColU","mamPriV","palAntN","proCap-Pcap_2.0_HiC","triManLat2","choHof-C_hoffmanni-2.0.1_HiC","chrAsi1m","dasNov3","eleEdw1m", "oryAfe2", "echTel2"]
    return expand(search_path, genome = RBBGenomes)
    
    
def getGeneCopyTableRBB_Atlantogenata_hq(wildcards):
    search_path = "output/geneCopyTable/{searchType}".format(searchType=wildcards.searchType) + "/{genome}_geneCopyTable_RBB_filtered.tsv"
    RBBGenomes = ["eleMaxD","loxAfr4","loxCycF","proCap-Pcap_2.0_HiC","triManLat2","choHof-C_hoffmanni-2.0.1_HiC","dasNov3", "oryAfe2", "echTel2"]
    return expand(search_path, genome = RBBGenomes)
    
    
def getGeneCopyTableEvidenced(wildcards):
    try:
        search_path = "output/geneCopyTable/{searchType}".format(searchType=wildcards.searchType) + "/{genome}_geneCopyTable_evidenced_filtered.tsv"
        RBBGenomes, = glob_wildcards(search_path)
        files = expand(search_path, genome = RBBGenomes)
        if files:
            return files
        else:
            raise Exception("No files found! Files: {0}\nSearch Path: {1}".format(str(files), search_path))
    except Exception:
        raise Exception("No files found! Files: {0}\nSearch Path: {1}".format(str(files), search_path))
    
rule consolidateGeneCopyTableRBB:
    input: 
        geneCopyTables = getGeneCopyTableRBB,
        genomeTable = "output/other/genomeTable.csv"
    output: 
        longTable = "output/geneCopyTable/{searchType}/consolidated-GeneCopyTable_RBB_filtered_long.csv",
        wideTable = "output/geneCopyTable/{searchType}/consolidated-GeneCopyTable_RBB_filtered_wide.csv",
        codedphylip = "output/phylip/{searchType}/consolidated-GeneCopyTable_RBB_filtered_coded.phylip",
        phylipGeneList = "output/phylip/{searchType}/consolidated-GeneCopyTable_RBB_filtered_coded_list.txt",
        wideTable_dyn = "output/geneCopyTable/{searchType}/consolidated-GeneCopyTable_RBB_filtered_wide_dyn.csv",
        codedphylip_dyn = "output/phylip/{searchType}/consolidated-GeneCopyTable_RBB_filtered_coded_dyn.phylip"
    script: "../code/consolidateGeneCopyTable.R"

rule consolidateGeneCopyTableEvidenced:
    input: 
        geneCopyTables = getGeneCopyTableEvidenced,
        genomeTable = "output/other/genomeTable.csv"
    output: 
        longTable = "output/geneCopyTable/{searchType}/eutheria-GeneCopyTable_evidenced_filtered_long.csv",
        wideTable = "output/geneCopyTable/{searchType}/eutheria-GeneCopyTable_evidenced_filtered_wide.csv",
        wideTable_dyn = "output/geneCopyTable/{searchType}/consolidated-GeneCopyTable_evidenced_filtered_wide_dyn.csv",
        codedphylip_dyn = "output/phylip/{searchType}/consolidated-GeneCopyTable_evidenced_filtered_coded_dyn.phylip",
        codedphylip = "output/phylip/{searchType}/consolidated-GeneCopyTable_evidenced_filtered_coded.phylip",
        phylipGeneList = "output/phylip/{searchType}/consolidated-GeneCopyTable_evidenced_filtered_coded_list.txt"
    script: "../code/consolidateGeneCopyTable.R"

rule consolidateGeneCopyTableRBB_Atlantogenata:
    input: 
        geneCopyTables = getGeneCopyTableRBB_Atlantogenata,
        genomeTable = "output/other/genomeTable.csv"
    output: 
        longTable = "output/geneCopyTable/{searchType}/atlantogenata-GeneCopyTable_RBB_filtered_long.csv",
        wideTable = "output/geneCopyTable/{searchType}/atlantogenata-GeneCopyTable_RBB_filtered_wide.csv",
        wideTable_dyn = "output/geneCopyTable/{searchType}/atlantogenata-GeneCopyTable_RBB_filtered_wide_dyn.csv",
        codedphylip = "output/phylip/{searchType}/atlantogenata-GeneCopyTable_RBB_filtered_coded.phylip",
        codedphylip_dyn = "output/phylip/{searchType}/atlantogenata-GeneCopyTable_RBB_filtered_coded_dyn.phylip",
        phylipGeneList = "output/phylip/{searchType}/atlantogenata-GeneCopyTable_RBB_filtered_coded_list.txt"
    script: "../code/consolidateGeneCopyTable.R"
    
rule consolidateGeneCopyTableEvidenced_Atlantogenata:
    input: 
        geneCopyTables = getGeneCopyTableEvidenced_Atlantogenata,
        genomeTable = "output/other/genomeTable.csv"
    output: 
        longTable = "output/geneCopyTable/{searchType}/atlantogenata-GeneCopyTable_evidenced_filtered_long.csv",
        wideTable = "output/geneCopyTable/{searchType}/atlantogenata-GeneCopyTable_evidenced_filtered_wide.csv",
        wideTable_dyn = "output/geneCopyTable/{searchType}/atlantogenata-GeneCopyTable_evidenced_filtered_wide_dyn.csv",
        codedphylip = "output/phylip/{searchType}/atlantogenata-GeneCopyTable_evidenced_filtered_coded.phylip",
        codedphylip_dyn = "output/phylip/{searchType}/atlantogenata-GeneCopyTable_evidenced_filtered_coded_dyn.phylip",
        phylipGeneList = "output/phylip/{searchType}/atlantogenata-GeneCopyTable_evidenced_filtered_coded_list.txt"
    script: "../code/consolidateGeneCopyTable.R"
    

# rule consolidateGeneCopyTableRBB_Atlantogenata_HQ:
#     input: 
#         geneCopyTables = getGeneCopyTableRBB_Atlantogenata_hq,
#         genomeTable = "output/other/genomeTable.csv"
#     output: 
#         longTable = "output/geneCopyTable/{searchType}/atlantogenata_hq-GeneCopyTable_RBB_filtered_long.csv",
#         wideTable = "output/geneCopyTable/{searchType}/atlantogenata_hq-GeneCopyTable_RBB_filtered_wide.csv",
#         wideTable_dyn = "output/geneCopyTable/{searchType}/atlantogenata_hq-GeneCopyTable_RBB_filtered_wide_dyn.csv",
#         codedphylip = "output/phylip/{searchType}/atlantogenata_hq-GeneCopyTable_RBB_filtered_coded.phylip",
#         codedphylip_dyn = "output/phylip/{searchType}/atlantogenata_hq-GeneCopyTable_RBB_filtered_coded_dyn.phylip",
#         phylipGeneList = "output/phylip/{searchType}/atlantogenata_hq-GeneCopyTable_RBB_filtered_coded_list.txt"
#     script: "../code/consolidateGeneCopyTable.R"
    
rule portTableToPublicationTable:
    input: 
        portTable="data/portTable.csv"
    output: 
        genomeTable = "output/other/genomeTable.csv"
    conda: "../envs/conda_R2.yaml"
    script: "../code/portTableToPublicationTable.R"
    
rule parsimony:
    input:
        tree = "data/stableTraits/{clade}.tree",
        wideTable = "output/geneCopyTable/{searchType}/{clade}-GeneCopyTable_RBB_filtered_wide{dyn}.csv"
    output: "output/parsimony/{searchType}/{clade}_RBB_filtered-parsimony{dyn,(_dyn)?}.nexus"
    conda: "../envs/conda_R2.yaml"
    script: "../code/parsimony.R"


rule analyzeParsimony:
    input: "output/parsimony/{searchType}/{clade}_RBB_filtered-parsimony.nexus"
    output: directory("output/geneLists/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/{clade}_parsimony_filtered/")
    conda: "../envs/conda_R2.yaml"
    script: "../code/analyzeParsimony.R"
    
rule maximumLikelihood: 
    input: 
        phylip = "output/phylip/{searchType}/{clade}-GeneCopyTable_{type}_filtered_coded{dyn}.phylip",
        #output/phylip/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata-GeneCopyTable_RBB_filtered_coded.phylip
        tree = "data/stableTraits/{clade}.tree" 
        #data/stableTraits/atlantogenata.tree
    output:
        log = "output/iqtree/{searchType}/maxLikelihood_model_{model}-dataType_{st}-asrMin_{asr_min}/{clade}_{type,RBB|evidenced}_filtered{dyn,(_dyn)?}.log",
        checkpoint = "output/iqtree/{searchType}/maxLikelihood_model_{model}-dataType_{st}-asrMin_{asr_min}/{clade}_{type,RBB|evidenced}_filtered{dyn,(_dyn)?}.ckp.gz",
        iqtree = "output/iqtree/{searchType}/maxLikelihood_model_{model}-dataType_{st}-asrMin_{asr_min}/{clade}_{type,RBB|evidenced}_filtered{dyn,(_dyn)?}.iqtree",
        treefile = "output/iqtree/{searchType}/maxLikelihood_model_{model}-dataType_{st}-asrMin_{asr_min}/{clade}_{type,RBB|evidenced}_filtered{dyn,(_dyn)?}.treefile",
        states = "output/iqtree/{searchType}/maxLikelihood_model_{model}-dataType_{st}-asrMin_{asr_min}/{clade}_{type,RBB|evidenced}_filtered{dyn,(_dyn)?}.state"
    params:
        outputDir = "output/iqtree/{searchType}/maxLikelihood_model_{model}-dataType_{st}-asrMin_{asr_min}/",
        prefix = "output/iqtree/{searchType}/maxLikelihood_model_{model}-dataType_{st}-asrMin_{asr_min}/{clade}_{type}_filtered{dyn,(_dyn)?}"
    threads: 24
    # shadow: "shallow"
    conda: "../envs/conda_iqtree.yaml"
    #shell: "iqtree -s Atlantogenata.phy -st MORPH -nt AUTO -m MK+FQ+I+G4 -te Atlantogenata.tree -asr -asr-min 0.8"
    shell: "mkdir -p {params.outputDir}; iqtree -s {input.phylip} -te {input.tree} -nt {threads} -st {wildcards.st} -m {wildcards.model} -asr -asr-min {wildcards.asr_min} -pre {params.prefix}"
    

rule analyzeMaximumLikelihood_dyn: 
    input: 
        iqtree = "output/iqtree/{searchType}/maxLikelihood_{modelType}/{clade}_{type}_filtered_dyn.treefile",
        states = "output/iqtree/{searchType}/maxLikelihood_{modelType}/{clade}_{type}_filtered_dyn.state",
        wideTable_dyn = "output/geneCopyTable/{searchType}/{clade}-GeneCopyTable_{type}_filtered_wide_dyn.csv",
        sitenames = "output/phylip/{searchType}/atlantogenata_hq-GeneCopyTable_RBB_filtered_coded_list.txt"
        #data/stableTraits/atlantogenata.tree
    output: 
        nexus = "output/geneLists/maxLikelihood_{modelType}/{searchType}/{clade}_{type,RBB|evidenced}_filtered_dyn/{clade}_{type}_filtered_dyn.nexus"
    params:
        dir = "output/geneLists/maxLikelihood_{modelType}/{searchType}/{clade}_{type}_filtered_dyn/"
    threads: 24
    conda: "../envs/conda_R2.yaml"
    script: "../code/analyzeML_dyn.R"
    
    
rule RBHBTranscriptTable:
    input: 
        rbhb_bed = "output/finalBed/{searchType}/{genome}_{type}_final.bed",
        transcripts_bed = "output/StringTie/finalMerge/{genome}-finalGuide.gff3"
    output: "output/RBHBID-to-StringTieTranscriptID/{searchType}/{genome}_{type}_ID-TranscriptID.tsv"
    conda: "../envs/conda_blat.yaml"
    shell:
        "tail -n+2 {input.rbhb_bed} | bedtools intersect -a stdin -b {input.transcripts_bed} -wa -wb | cut -f 4,21 | grep -v exon_number > {output}"
