rule prep_TOGA_gff:
    input: 
        gff = '../data/toga/{genome}/{ref}/query_annotation.agat.gff',
        geneids = '../data/toga/{genome}/{ref}/query_isoforms.tsv'
    output: '../data/toga/{genome}/{ref}/query_annotation.agat.asMRNA.gff3'
    script: "../code/prep_AGAT_GFF.R"

rule TOGA_genespan_bed2gff:
    input: 
        bed = '../data/toga/{genome}/{ref}/query_gene_spans.bed'
    output:
        qsgff = '../data/toga/{genome}/{ref}/query_gene_spans.gff3',
    shell: "agat_convert_bed2gff.pl --bed {input.bed} -o {output.qsgff}"

rule TOGA_genespan_gff_fixnames:
    input: 
        qsgff = '../data/toga/{genome}/{ref}/query_gene_spans.gff3',
    output:
        cngff = '../data/toga/{genome}/{ref}/query_gene_spans.correctedName.gff3',
    shell: 
        " cat {input.qsgff} | "
        "  sed 's/\\tdata\\t/\\tTOGA\\t/; s/ID=[0-9]\+;Name=\([A-Za-z_0-9]\+\)/ID=\\1;Name=gene.\\1/'"
        "  > {output.cngff}"

rule TOGA_regGene_clean:
    input: 
        cngff = '../data/toga/{genome}/{ref}/query_gene_spans.correctedName.gff3',
        gff = '../data/toga/{genome}/{ref}/query_annotation.agat.asMRNA.gff3'
    output:
        reg_gff = '../data/toga/{genome}/{ref}/query_annotation.agat.asMRNA.mergedwithreg_gene.gff3',
        reg_clean_gff = '../data/toga/{genome}/{ref}/query_annotation.agat.asMRNA.mergedwithreg_gene.clean.gff3'
    shell: 
        " cat {input.cngff} {input.gff} > {output.reg_gff}; "
        " agat_convert_sp_gxf2gxf.pl -g {output.reg_gff} -o {output.reg_clean_gff}"


rule TOGA_RegsOnly:
    input:
        reg_clean_gff = '../data/toga/{genome}/{ref}/query_annotation.agat.asMRNA.mergedwithreg_gene.clean.gff3'
    output:
        regs = '../data/toga/{genome}/{ref}/query_annotation.agat.asMRNA.mergedwithreg_gene.clean.gff3.regs',
        reg_clean_gff = '../data/toga/{genome}/{ref}/query_annotation.regsOnly.gff3'
    shell:
        "grep 'reg_' {input.reg_clean_gff} | "
        " grep -oPe '(?<=ID=)reg_[0-9]+(?=;)' > {output.regs}; "
        " agat_sp_filter_feature_from_keep_list.pl "
        " --gff {input.reg_clean_gff} "
        " --keep_list {output.regs} "
        " -o {output.reg_clean_gff}"

rule TOGA_APPRISOnly:
    input:
        appris = 'data/appris.principalOnly.ENST.txt',
        reg_clean_gff = '../data/toga/{genome}/{ref}/query_annotation.regsOnly.gff3'
    output:
        mrna = '../data/toga/{genome}/{ref}/query_annotation.regsOnly.clean.mRNAs',
        genes = '../data/toga/{genome}/{ref}/query_annotation.regsOnly.clean.APPRIS.txt',
        gff = '../data/toga/{genome}/{ref}/query_annotation.regsOnly.clean.APPRIS.gff3'
    threads: 32
    shell:
        "grep -Pe '\\tmRNA\\t' {input.reg_clean_gff} > {output.mrna}; "
        " cat {input.appris} | "
        " parallel -j {threads} "
        " 'grep {{}} {output.mrna}' | "
        " grep -oPe '(?<=Parent=)[-._A-Za-z0-9]+(?=;)' | "
        " sort -u > {output.genes}; "
        " agat_sp_filter_feature_from_keep_list.pl "
        " --gff {input.reg_clean_gff} "
        " --keep_list {output.genes} "
        " -o {output.gff}"

rule TOGA_filter_incomplete:
    input: 
        gff = '../data/toga/{genome}/{ref}/{gff}.gff3',
        genome = 'data/genomes/{genome}.fa'
    output: '../data/toga/{genome}/{ref}/{gff}.completeCDS.gff3'
    shell: 
        "agat_sp_filter_incomplete_gene_coding_models.pl "
        " --gff {input.gff} "
        " --fasta {input.genome} "
        " -o {output}"

rule EVM_gather_evidence:
    input:
        augustus = 'output/funannotate/{genome}_TOGA_Prot/predict_misc/augustus.evm.gff3',
        liftoff = 'data/liftoff/mMyoMyo1-{genome}-GenBank.mRNA.withCDS.cdsphase.flagPseudogenes.gff_polished',
        toga = '../data/toga/{genome}/hg38/query_annotation.regsOnly.gff3',
        miniprot = 'output/miniprot_0.6-r194-dirty/{genome}_uniparc-myotis-chiroptera.fixedinheritance.merged.gff',
        genome_softmasked = 'output/funannotate/{genome}_TOGA_Prot/predict_misc/genome.softmasked.fa',
        repeatmasker = '../data/repeatmasker/union/{genome}.full.denovo.gff',
        protein_alignment = 'output/funannotate/{genome}_TOGA_Prot/predict_misc/protein_alignments.gff3',
        transcript_alignment = 'output/funannotate/{genome}_TOGA_Prot/predict_misc/transcript_alignments.gff3',
        weights = 'data/weights.evm.txt'
    output:
        augustus = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_augustus.gff3',
        liftoff = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_liftoff.gff3',
        toga = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_TOGA.gff3',
        miniprot = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_miniprot.gff3',
        genome_softmasked = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/genome.softmasked.fa',
        repeatmasker = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/repeatmasker.gff',
        protein_alignment = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/protein_alignments.gff3',
        transcript_alignment = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/transcript_alignments.gff3',
        weights = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/weights.evm.txt'
    run:
        shell('mkdir output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/'.format(genome=genome))
        for i,j in zip(input,output):
            shell('cp {} {}'.format(i,j))

rule EVM_merge_evidence:
    input:
        augustus = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_augustus.gff3',
        liftoff = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_liftoff.gff3',
        toga = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_TOGA.gff3',
        miniprot = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/evidence_miniprot.gff3'
    output: 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/gene_predictions.gff3'
    shell: 'cat {input} > {output}'

rule EVM_run:
    input:
        genes = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/gene_predictions.gff3',
        genome_softmasked = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/genome.softmasked.fa',
        repeatmasker = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/repeatmasker.gff',
        protein_alignment = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/protein_alignments.gff3',
        transcript_alignment = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/transcript_alignments.gff3',
        weights = 'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/weights.evm.txt'
    output: 
        'output/EVM_manual/{genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/{genome}_TOGARegOnly_Augustus_MiniprotMerged_GenBankLiftoff.EVM.gff3'
    threads: 40
    shell:
        'cd output/EVM_manual/{wildcards.genome}_TOGARegOnly_Augustus_MiniprotMergedCompleteCDS_GenBankLiftoff/; '
        'singularity exec ../../../evidencemodeler_latest.sif '
        ' EVidenceModeler --debug --CPU {threads} '
        '--segmentSize 1000000 '
        '--overlapSize 150000 '
        '--sample_id {wildcards.genome}_TOGARegOnly_Augustus_MiniprotMerged_GenBankLiftoff '
        '--genome genome.softmasked.fa '
        '--weights weights.evm.txt '
        '--gene_predictions gene_predictions.gff3 '
        '--protein_alignments protein_alignments.gff3 '
        '--repeats repeatmasker.gff'


rule EVM_diamond_keeper:
    input: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond.tsv'
    output:
        tbl = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond.criteria.tsv',
        keep = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond.keep.ids'
    script: '../code/diamond_keepers.R'


rule EVM_diamond_pass:
    input: 
        gff = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.gff3',
        keep = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond.keep.ids'        
    output: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond_pass.gff3'
    shell: 'agat_sp_filter_feature_from_keep_list.pl --gff {input.gff} --keep_list {input.keep} -o {output}'

rule EVM_diamond_fail:
    input: 
        gff = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.gff3',
        keep = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond.keep.ids'        
    output: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond_fail.gff3'
    shell: 'agat_sp_filter_feature_from_kill_list.pl --gff {input.gff} --kill_list {input.keep} -o {output}'

rule EVM_diamond_fail_keep_BigORF:
    input: 
        gff = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond_fail.gff3',
    output: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond_fail.CDSOver120_sup=120.gff'
    params:
        orf = 120,
        prefix = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond_fail.CDSOver120.gff'
    shell: 'agat_sp_filter_by_ORF_size.pl --gff {input.gff} --test ">=" -s {params.orf} -o {params.prefix}'
    
rule EVM_diamond_final_pass:
    input: 
        gff1 = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond_pass.gff3',
        gff2 = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.diamond_fail.CDSOver120_sup=120.gff'
    output: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.final_pass.gff3'
    shell: 'agat_sp_merge_annotations.pl -f {input.gff1} -f {input.gff2} --out {output}'


rule EVM_finalpass_getCDS:
    input: 
        gff = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.final_pass.gff3',
        fasta = 'data/genomes/{genome}.fa'
    output: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.final_pass.cds.fa'
    wildcard_constraints:
        genome = '[A-Za-z]+[0-9]'
    shell: 'agat_sp_extract_sequences.pl --gff {input.gff} -t cds -f {input.fasta} --out {output}'


rule EVM_finalpass_getProtein:
    input: 
        gff = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.final_pass.gff3',
        fasta = 'data/genomes/{genome}.fa'
    output: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.final_pass.faa'
    wildcard_constraints:
        genome = '[A-Za-z]+[0-9]'
    shell: 'agat_sp_extract_sequences.pl --gff {input.gff} -p -f {input.fasta} --out {output}'

rule EVM_finalpass_flagPseudogenes:
    input: 
        gff = 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.final_pass.gff3',
        fasta = 'data/genomes/{genome}.fa'
    output: 'output/EVM_manual/{genome}_{evidence}/{genome}_{evm}.complementedTOGA.final_pass.flagPseudogenes.gff3'
    wildcard_constraints:
        genome = '[A-Za-z]+[0-9]'
    shell: 'agat_sp_flag_premature_stop_codons.pl -f {input.fasta} --gff {input.gff} -o {output}'
        
