rule miniprot_index:
    input: "data/genomes/{genome}.fa"
    output: "data/genomes/{genome}.mpi"
    threads: 32
    shell: "code/miniprot/miniprot -t{threads} -d {output} {input}"

rule miniprot_gff:
    version: "0.6-r194-dirty"
    input: 
        idx = "data/genomes/{genome}.mpi",
        prot = "data/protein_evidence/{proteins}.fa"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.gff"
    threads: 32
    shell: " if [ $(code/miniprot/miniprot --version) != '0.6-r194-dirty' ]; then print 'wrong version'; else code/miniprot/miniprot -ut{threads} --gff {input.idx} {input.prot} > {output}; fi"

rule miniprot_cdsfix_gff:
    version: "0.6-r194-dirty"
    input: 
        gff = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.gff",
        fa = "data/genomes/{genome}.fa"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.cdsfix.gff"
    threads: 32
    shell: " if [ $(code/miniprot/miniprot --version) != '0.6-r194-dirty' ]; "
           " then print 'wrong version'; else "
           " ~/miniconda3/envs/AGAT/bin/perl ~/miniconda3/envs/AGAT/bin/agat_sp_fix_cds_phases.pl --gff {input.gff} --fasta {input.fa} -o {output}; "
           " fi"

rule miniprot_fix_inheritance:
    input: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.gff"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.gff"
    script: "../code/fix_inheritance_gff.py"

rule miniprot_getFASTA:
    input: 
        gff="output/miniprot_0.6-r194-dirty/{genome}_{variant}.gff",
        genome="data/genomes/{genome}.fa"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{variant}.faa"
    shell: 'agat_sp_extract_sequences.pl -f {input.genome} -p -g {input.gff} -o {output} '

rule miniprot_gtf:
    version: "0.6-r194-dirty"
    input: 
        idx = "data/genomes/{genome}.mpi",
        prot = "data/protein_evidence/{proteins}.fa"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.gtf"
    threads: 32
    shell: " if [ $(code/miniprot/miniprot --version) != '0.6-r194-dirty' ]; then print 'wrong version'; else code/miniprot/miniprot -ut{threads} --gtf {input.idx} {input.prot} > {output}; fi"

    
rule miniprot_format_EVM:
    """
    EVM wants a parent-less GFF3 for protein evidence where parentage is implied based on
    sharing the ID field.
    """
    version: 2
    input: "output/miniprot_{version}/{genome}_{proteins}.gff"
    output: "output/miniprot_{version}/{genome}_{proteins}.evm.gff"
    shell: 
        "cat {input} | "
        " grep CDS | "
        " sed 's/Parent=/ID=/; s/CDS/protein_to_nucleotide_match/' > "
        " {output}"


rule miniprot_collapse_gff:
    version: "0.6-r194-dirty"
    input: 
        gff="output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.gff",
        fai = "data/genomes/{genome}.fa.fai"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.merged.gff"
    threads: 32
    params:
           tmp = "output/miniprot_0.6-r194-dirty/temp_{genome}_{proteins}",
           outprefix = "output/miniprot_0.6-r194-dirty/temp_{genome}_{proteins}/{genome}_{proteins}.fixedinheritance"
    shell: " if [ $(code/miniprot/miniprot --version) != '0.6-r194-dirty' ]; "
           " then print 'wrong version'; else "
           " mkdir {params.tmp} && echo 'Created {params.tmp}'; "
           " cut -f1 {input.fai} | "
           "   parallel 'grep {{}} {input.gff} > {params.outprefix}.{{}}.gff' "
           " && echo 'Splitted GFF by scaffold'; "
           " find {params.outprefix}.* -type f -empty -delete "
           " && echo 'Removed empty scaffold GFFs'; "
           " ls {params.outprefix}.*.gff | "
           "   parallel 'agat_sp_fix_overlaping_genes.pl "
           "     --gff {{}} -o {{.}}.merged.gff' "
           " && echo 'Collapsed genes in each scaffold'; "
           " cat {params.outprefix}.*.merged.gff > {output} "
           " && echo 'Merged into final GFF'; "
           " rm -r {params.tmp} && echo 'Cleaned up.'; "
           " fi"

rule miniprot_completeCDS_gff:
    version: "0.6-r194-dirty"
    input: 
        gff = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.merged.gff",
        fa = "data/genomes/{genome}.fa",
        fai = "data/genomes/{genome}.fa.fai"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.merged.completeCDS.gff"
    threads: 32
    params:
           tmp = "output/miniprot_0.6-r194-dirty/temp_{genome}_{proteins}_completeCDS",
           outprefix = "output/miniprot_0.6-r194-dirty/temp_{genome}_{proteins}_completeCDS/{genome}_{proteins}.fixedinheritance.merged"
    shell: " if [ $(code/miniprot/miniprot --version) != '0.6-r194-dirty' ]; "
           " then print 'wrong version'; else "
           " mkdir {params.tmp}; "
           " cut -f1 {input.fai} | "
           "   parallel 'grep {{}} {input.gff} > {params.outprefix}.{{}}.gff'; "
           " find {params.outprefix}.* -type f -empty -print -delete; "
           " ls {params.outprefix}.*.gff | "
           "   parallel 'agat_sp_filter_incomplete_gene_coding_models.pl "
           "     --gff {{}} --fasta {input.fa} -o {{.}}.completeCDS.gff'; "
           " cat {params.outprefix}.*.completeCDS.gff > {output}; "
           " rm -r {params.tmp}; "
           " fi"

rule miniprot_completeCDS_clean:
    input: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.merged.completeCDS.gff"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.merged.completeCDS.clean.gff"
    shell: "agat_convert_sp_gxf2gxf.pl -g {input} -o {output}"

#rule miniprot_RBH_getFAA:
#    input: 
#        gff = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.gff",
#        genome = "data/genomes/{genome}.fa"
#    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.faa"
#    shell: 'agat_sp_extract_sequences.pl -f {input.genome} -p -g {input.gff} -o {output} '

rule miniprot_RBH_DIAMOND:
    input: 
        faa = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.faa",
        db = "data/protein_evidence/{proteins}.dmnd"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.diamond.tsv"
    threads: 32
    shell: 
        "code/diamond blastp "
        " --query {input.faa} "
        " --db {input.db} "
        " --outfmt 6 qseqid bitscore sseqid pident length mismatch gapopen qlen qstart qend slen sstart send ppos evalue "
        " --ultra-sensitive --max-target-seqs 1 --evalue 1e-10 "
        " --threads {threads} "
        " > {output}"

rule miniprot_RBH_getHeader:
    input: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.faa"
    output: "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.header"
    shell: "grep '>' {input} | sed 's/>//' > {output}"

rule miniprot_RBH_getRBH:
    input: 
        diamond = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.diamond.tsv",
        gff =  "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.gff"
    output: 
        tsv = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.diamond_rbh.tsv",
        txt = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.diamond_rbh.keep"
    script: '../code/diamond_getRBH.R'
    
rule miniprot_RBH_getGFF:
    input:
        gff =  "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.gff",
        txt = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.diamond_rbh.keep"
    output:
        gff = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.diamond_rbh.gff3",
        report = "output/miniprot_0.6-r194-dirty/{genome}_{proteins}.fixedinheritance.diamond_rbh_report.txt"
    shell: 
        "if [ -f {output.report} ]; then rm {output.report}; fi; "
        "if [ -f {output.gff} ]; then rm {output.gff}; fi; "
        "agat_sp_filter_feature_from_keep_list.pl "
        "--gff {input.gff} --keep_list {input.txt} -o {output.gff}"
