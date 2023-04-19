rule EVM_drop_notmRNA:
    input: 'data/liftoff/{ref}-{genome}-{target}.gff_polished'
    output: 
        notmRNA = 'data/liftoff/{ref}-{genome}-{target}.gff_polished.notmRNA.txt',
        gff='data/liftoff/{ref}-{genome}-{target}.mRNA.gff_polished'
    threads: 40
    shell: "grep -vPe '\\t(mRNA|exon|gene|CDS)\\t' {input} | "
           " grep -oPe '(?<=ID=)[-._A-Za-z0-9]+(?=;)' | "
           "sort -u > {output.notmRNA}; "
           " if [ -z '$(cat {output.notmRNA})' ]; then "
           "cat {input} > {output.gff}; else "
           " cat {input} | "
           "  parallel -j{threads} --cat 'grep -vf {output.notmRNA} {{}}' > {output.gff}; "
           " fi; exit 0"

rule EVM_find_mRNA_missingCDS:
    input: 'data/liftoff/{ref}-{genome}-{target}.mRNA.gff_polished'
    output: 'data/liftoff/{ref}-{genome}-{target}.mRNA.nocds.txt'
    script: '../code/EVM_check_missing_cds.py'

rule EVM_drop_mRNA_missingCDS:
    input: 
        gff='data/liftoff/{ref}-{genome}-{target}.mRNA.gff_polished',
        txt='data/liftoff/{ref}-{genome}-{target}.mRNA.nocds.txt'
    output: 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.gff_polished'
    threads: 40
    shell: 
        "if [ -z '$(cat {input.txt})' ]; then "
        "cat {input.gff} > {output}; else "
        "cat {input.gff} | parallel -j{threads} --cat 'grep -vf {input.txt} {{}}' > {output};"
        " fi"
        
rule agat_fixCDS:
    input: 
        gff='data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.gff_polished',
        genome = 'data/genomes/{genome}.fa'
    output: 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.gff_polished'
    shell: 'agat_sp_fix_cds_phases.pl -g {input.gff} -fa {input.genome} -o {output}'

rule agat_flagPseudogenes:
    input: 
        gff='data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.gff_polished',
        genome = 'data/genomes/{genome}.fa'
    output: 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.flagPseudogenes.gff_polished'
    shell: 'agat_sp_flag_premature_stop_codons.pl --gff {input.gff} --fasta {input.genome} --out {output}'

rule EVM_find_pseudogenes:
    input: 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.flagPseudogenes.gff_polished',
    output: 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.Pseudogenes.txt'
    shell: "grep ';pseudo=' {input} | grep -Pe '\\t(gene|mRNA)\\t' | grep -oPe '(?<=ID=)[A-Za-z]+[-_.A-Za-z0-9]+;' > {output}"

rule EVM_drop_pseudogenes:
    input: 
        gff='data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.flagPseudogenes.gff_polished',
        txt='data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.Pseudogenes.txt'
    output: 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.noPseudogenes.gff_polished'
    threads: 40
    shell: "cat {input.gff} | parallel -j{threads} --cat 'grep -vf {input.txt} {{}}' > {output}"

rule liftoff_getFasta:
    input: 
        genome = 'data/genomes/{genome}.fa',
        gff = 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.flagPseudogenes.gff_polished'
    output: 'data/liftoff/{ref}-{genome}-{target}.mRNA.withCDS.cdsphase.flagPseudogenes.faa'
    shell: 'agat_sp_extract_sequences.pl -f {input.genome} -p -g {input.gff} -o {output}'
