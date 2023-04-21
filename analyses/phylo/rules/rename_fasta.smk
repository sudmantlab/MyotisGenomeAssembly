GENOMES = ["mMyoAui1","mMyoLuc1","mMyoEvo1","mMyoThy1","mMyoVol1","mMyoVel1","mMyoCai1","mMyoOcc1","mMyoYum1"]


def get_BUSCO_union_ids():
    with open('data/genes/BUSCO_eutheria/busco_table/union_complete.list') as infile:
        return [i.strip() for i in infile]

def get_BUSCO_regroup_filenames_codon(wildcards):
    filepath = 'data/genes/BUSCO_eutheria/codon_nuc/byBUSCO/{sequence}.fa'
    buscos = get_BUSCO_union_ids()
    return [filepath.format(sequence = b) for b in buscos]

def get_BUSCO_regroup_filenames_prot(wildcards):
    filepath = 'data/genes/BUSCO_eutheria/prot_aa/byBUSCO/{sequence}.fa'
    buscos = get_BUSCO_union_ids()
    return [filepath.format(sequence = b) for b in buscos]    

rule regroup_BUSCO_complete_sequence_all_prot:
    input: 
        list = 'data/genes/BUSCO_eutheria/busco_table/union_complete.list', 
        buscos = get_BUSCO_regroup_filenames_prot

rule regroup_BUSCO_complete_sequence_all_codon:
    input: 
        list = 'data/genes/BUSCO_eutheria/busco_table/union_complete.list', 
        buscos = get_BUSCO_regroup_filenames_codon


rule get_BUSCO_complete_list:
    input: 
        'data/genomes/BUSCO_eutheria/busco_table/{genome}_complete.tsv'
    output:
        'data/genes/BUSCO_eutheria/busco_table/{genome}_complete.list'
    shell: 
        "cut -f2 {input} | sed 's/^/>/; s/$/ /' | sort -u > {output}"

rule union_BUSCO_complete_sequence_ids:
    input:
        list = expand('data/genes/BUSCO_eutheria/busco_table/{genome}_complete.tsv', genome=GENOMES)
    output:
        list = 'data/genes/BUSCO_eutheria/busco_table/union_complete.list'    
    run:
        sets = []
        for l in input.list:
            with open(l) as infile:
                sets.append([i.strip().split('\t')[0] for i in infile])
        union_BUSCO = list(set().union(*sets))
        with open(output.list, 'w') as outfile:
            outfile.write('\n'.join(union_BUSCO))

### REQUIRES 1-LINE FASTAs
### EG:
### ```
### cat mMyoCai1_TOGARegOnly_Augustus_MiniprotMerged_GenBankLiftoff.complementedTOGA.final_pass.faa | 
###   awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > 
###   data/genes/all/prot_aa/mMyoCai1.faa
### ```
    
rule get_BUSCO_complete_sequence:
    input: 
        fasta = 'data/genes/all/{type}/{genome}.fa',
        list = 'data/genes/BUSCO_eutheria/busco_table/{genome}_complete.list'
    output:
        'data/genes/BUSCO_eutheria/{type}/{genome}_complete.fa'
    shell: 
        'grep -A1 -f {input.list} {input.fasta} > {output}'

rule rename_fasta_BUSCO:
    input: 
        fasta = 'data/genes/BUSCO_eutheria/{type}/{genome}_complete.fa',
        tt = 'data/genes/BUSCO_eutheria/busco_table/{genome}_complete.tsv'
    output:
        'data/genes/BUSCO_eutheria/{type}/{genome}_complete.renamed.fa'
    script: 
        '../code/rename_fasta.py'


rule regroup_BUSCO_complete_sequence:
    input: 
        fasta = expand('data/genes/BUSCO_eutheria/{{type}}/{genome}_complete.renamed.fa', genome = GENOMES),
    output:
        fasta = 'data/genes/BUSCO_eutheria/{type}/byBUSCO/{sequence}.fa'
    shell:
        #'if [ ! -d data/genes/BUSCO_eutheria/prot_aa/byBUSCO/ ]; '
        #' then mkdir data/genes/BUSCO_eutheria/prot_aa/byBUSCO/; fi; '
        'cat {input.fasta} | grep -A1 "_{wildcards.sequence} " > {output.fasta}'


