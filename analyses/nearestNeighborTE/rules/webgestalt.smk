rule WebGestalt_duplicates:
    input:
        geneSet = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_{type}_filtered.txt",
        refGeneSet = "output/recBlastDBPrep/{query}_geneList_filtered.txt"
    output:
        outputDirectory = directory("output/ORA/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/duplicates/{enrichmentDB}-FDR_{fdr}/{genome}_{type}_filtered")
    conda: "../envs/conda_R.yaml"
    script: "../code/runORA.R"
  
  
rule WebGestalt_intersection:
    input:
        geneSet = "output/geneLists/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/intersection-{intersection}_filtered.txt",
        refGeneSet = "output/recBlastDBPrep/{query}_geneList_filtered.txt"
    output:
        outputDirectory = directory("output/ORA/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/intersection/{enrichmentDB}-FDR_{fdr}/intersection-{intersection}_filtered")
    conda: "../envs/conda_R.yaml"
    script: "../code/runORA.R"
  
  
rule WebGestalt_oneShot:
    input:
        geneSet = "output/geneLists/oneShot/{list}_filtered.txt",
        refGeneSet = "output/recBlastDBPrep/AvA_geneList_filtered.txt"
    output: 
        outputDirectory = directory("output/ORA/oneShot/{enrichmentDB}-FDR_{fdr}/{list}_filtered")
    conda: "../envs/conda_R.yaml"
    script: "../code/runORA.R"
  
#   
# rule WebGestalt_subset:
#     input:
#         geneSet = "output/geneLists/{searchType}/{subType}/{list}.txt",
#         refGeneSet = "output/recBlastDBPrep/AvA_geneList_filtered.txt"
#     output: 
#         outputDirectory = directory("output/ORA/{subType}/{searchType}/{enrichmentDB}-FDR_{fdr}/{list}")
#     params:
#         rootDir= "output/ORA/{subType}/{searchType}/{enrichmentDB}-FDR_{fdr}"
#     conda: "../envs/conda_R.yaml"
#     script: "../code/runORA.R"
  
rule WebGestalt_ML:
    input:
        geneSet = "output/geneLists/maxLikelihood_{modelType}/{searchType}/{cladeType}/{list}.txt",
        refGeneSet = "output/recBlastDBPrep/AvA_geneList_filtered.txt"
    output: 
        outputDirectory = directory("output/ORA/maxLikelihood_{modelType}/{searchType}/{cladeType}/{enrichmentDB}-FDR_{fdr}/{list}")
    params:
        rootDir= "output/ORA/maxLikelihood_{modelType}/{searchType}/{cladeType}/{enrichmentDB}-FDR_{fdr}"
    conda: "../envs/conda_R.yaml"
    script: "../code/runORA.R"
