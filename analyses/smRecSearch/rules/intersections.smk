
rule GeneCopyTableRBB:
    input:
        TranslationTable = ancient("data/input/{query}.tsv"),
        RBB_file = ancient("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb"),
        ECNC_file = ancient("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/ECNC/{genome}_RecBlastOutput.bed.rbb.ecnc")
    output:
        copyTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_RBB.tsv",
        finalBed = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_RBB_final.bed",
        lociTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_lociTable_RBB.tsv",
        exonTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_exonTable_RBB.tsv",
        duplicateGeneList = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicateList/{genome}_geneDuplicates_RBB.txt"
    conda: "../envs/conda_R2.yaml"
    threads: 100
    log: "logs/GeneCopyTableRBB/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_RBB.log"
    script: "../code/getCopyTable.R"


rule GeneCopyTableRBB_filtered:
    input:
        geneListFiltered = ancient("output/recBlastDBPrep/{query}_geneList_filtered.txt"),
        copyTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_RBB.tsv",
        finalBed = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_RBB_final.bed",
        lociTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_lociTable_RBB.tsv",
        exonTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_exonTable_RBB.tsv",
        duplicateGeneList = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicateList/{genome}_geneDuplicates_RBB.txt"
    output:
        copyTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_RBB_filtered.tsv",
        lociTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_lociTable_RBB_filtered.tsv",
        finalBed = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_RBB_final_filtered.bed",
        exonTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_exonTable_RBB_filtered.tsv",
        duplicateGeneList = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicateList/{genome}_geneDuplicates_RBB_filtered.txt"
    conda: "../envs/conda_R2.yaml"
    threads: 100
    log: "logs/GeneCopyTableRBB_filtered/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}-{genome}_geneCopyTable_RBB_filtered.log"
    script: "../code/FilterCopyTable.R"

rule GeneCopyTableEvidenced:
    input:
        TranslationTable = ancient("data/input/{query}.tsv"),
        RBB_file = ancient("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_evidenced.bed"),
        ECNC_file = ancient("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/ECNC/{genome}_RecBlastOutput.bed.rbb.ecnc")
    output:
        copyTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_evidenced.tsv",
        lociTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_lociTable_evidenced.tsv",
        finalBed = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_evidenced_final.bed",
        exonTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_exonTable_evidenced.tsv",
        duplicateGeneList = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicateList/{genome}_geneDuplicates_evidenced.txt"
    conda: "../envs/conda_R2.yaml"
    threads: 100
    log: "logs/GeneCopyTableEvidenced/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_evidenced.log"
    script: "../code/getCopyTable.R"

rule GeneCopyTableEvidenced_filtered:
    input:
        geneListFiltered = ancient("output/recBlastDBPrep/{query}_geneList_filtered.txt"),
        copyTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_evidenced.tsv",
        lociTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_lociTable_evidenced.tsv",
        finalBed = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_evidenced_final.bed",
        exonTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_exonTable_evidenced.tsv",
        duplicateGeneList = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicateList/{genome}_geneDuplicates_evidenced.txt"
    output:
        copyTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_evidenced_filtered.tsv",
        lociTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_lociTable_evidenced_filtered.tsv",
        finalBed = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_evidenced_final_filtered.bed",
        exonTable = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_exonTable_evidenced_filtered.tsv",
        duplicateGeneList = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicateList/{genome}_geneDuplicates_evidenced_filtered.txt"
    conda: "../envs/conda_R2.yaml"
    threads: 100
    log: "logs/GeneCopyTableEvidenced_Filtered/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_evidenced_filtered.log"
    script: "../code/FilterCopyTable.R"

rule ReferenceGeneList:
    input:
        TranslationTable = "data/input/{query}.tsv"
    output:
        geneList = "output/recBlastDBPrep/{query}_geneList.txt",
        filteredGeneList = "output/recBlastDBPrep/{query}_geneList_filtered.txt"
    conda: "../envs/conda_R2.yaml"
    log: "logs/ReferenceGeneList/{query}_geneList.log"
    script: "../code/getReferenceGeneList.R"

rule linkGeneLists:
    input: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicateList/{genome}_geneDuplicates_{type}.txt"
    output: "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_geneDuplicates_{type}.txt"
    shell: "ln -sf ../../../{input} {output}"


rule linkCopyTable:
    input: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/geneCopyTable/{genome}_geneCopyTable_{type}.tsv"
    output: "output/geneCopyTable/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_geneCopyTable_{type}.tsv"
    shell: "ln -sf ../../../{input} {output}"

      
rule linkFinalBed:
    input: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_{type}.bed"
    output: "output/finalBed/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}.bed"
    shell: "ln -sf ../../../{input} {output}"


rule getDuplicatesBed:
    input:
        geneList = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_geneDuplicates_{type}_filtered.txt",
        bed = "output/finalBed/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_final_filtered.bed"
    output: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/finalBed/{genome}_{type}_duplicatesOnly_filtered.bed"
    shell: "grep -F -f {input.geneList} {input.bed} > {output}"

rule getIntersections:
    input:
        geneListA = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeA}_{typeA}_filtered.txt",
        geneListB = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genomeB}_{typeB}_filtered.txt"
    output:
        A_not_B = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/intersection-{genomeA}_{typeA}_filtered-not-{genomeB}_{typeB}_filtered.txt",
        B_not_A = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/intersection-{genomeB}_{typeB}_filtered-not-{genomeA}_{typeA}_filtered.txt",
        A_and_B = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/intersection-{genomeA}_{typeA}_filtered-and-{genomeB}_{typeB}_filtered.txt"
    conda: "../envs/conda_R.yaml"
    script: "../code/intersectGeneLists.R"

