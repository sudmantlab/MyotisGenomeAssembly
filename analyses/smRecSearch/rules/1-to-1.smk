rule get_1to1_hits:
    input: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb"
    output: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/1-to-1/{genome}_RecBlastOutput.bed.rbb"
    conda: "../envs/conda_R2.yaml"
    script: "../code/get_1-to-1.R"


