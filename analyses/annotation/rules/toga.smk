rule toga_to_gff_AGAT:
    input: '../data/toga/{genome}/{ref}/query_annotation.bed'
    output: '../data/toga/{genome}/{ref}/query_annotation.agat.gff'
    shell: 'agat_convert_bed2gff.pl --bed {input} --source TOGA --gff {output}'


rule toga_remove_UTRs:
    input: '../data/toga/{genome}/{ref}/query_annotation.agat.gff'
    output: '../data/toga/{genome}/{ref}/query_annotation.agat.noUTR.gff'
    shell: "cat {input} | sed '/prime_UTR\t/d' > {output}"


rule toga_fixFeatures:
    input: "../data/toga/{genome}/{ref}/query_annotation.agat.noUTR.gff"
    output: "../data/toga/{genome}/{ref}/query_annotation.agat.noUTR.mRNA.gff"
    #script: "../code/fixTOGAGFF.R"
    shell: "agat_convert_sp_gxf2gxf.pl -g {input} -o {output}"
 
rule liftoff_getFasta:
    input:
        genome = 'data/genomes/{genome}.fa',
        gff = "../data/toga/{genome}/{ref}/query_annotation.agat.noUTR.mRNA.gff"
    output: "../data/toga/{genome}/{ref}/query_annotation.agat.noUTR.mRNA.faa"
    shell: 'agat_sp_extract_sequences.pl -f {input.genome} -p -g {input.gff} -o {output}'
