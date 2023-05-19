rule ORATable:
    input:
        directory("output/ORA/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata_RBB_filtered_dyn/pathway_Reactome-top100/")
    output:
         "../output/pubFiles/ORA_table_pub.csv"
    # conda: "../envs/R-geneLists.yaml"
    script: "../code/pub_ORATable.R"


rule bodySizeTree:
    input:
        bsize = "output/AncSizeRec/AnnotatedTreeForFigure.tree",
        nexus = "data/stableTraits/eutheria.nexus",
        ancstates = "data/stableTraits/eutheria.ancstates",
        brlens = "data/stableTraits/eutheria.brlens",
        genomeTable = "output/other/genomeTable.csv"
    output:
        outTree_full = "output/pubFiles/bodysize-Eutheria.csv",
        outTree_Atlantogenata = "output/pubFiles/bodysize-Atlantogenata.csv"
    script: "../code/pub_bodySizeTree.R"


rule lifespanPGLS_cancerRisk:
    input:
        anage = "data/AnAge/anage_build14.txt",
        Atlantogenata = "output/pubFiles/bodysize-Atlantogenata.csv"
    output:
        cancerSuccep = "output/pubFiles/cancerRisk-Atlantogenata.csv",
        latex = "output/pubFiles/lifespan-PGLS.tex",
        html = "output/pubFiles/lifespan-PGLS.html"
    script: "../code/pub_lifespanPGLS_cancerRisk.R"


