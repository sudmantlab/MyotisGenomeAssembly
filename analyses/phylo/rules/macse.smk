def get_BUSCO_union_ids():
    with open('data/genes/BUSCO_eutheria/busco_table/union_complete.list') as infile:
        return [i.strip() for i in infile]

def get_BUSCO_union_ids_glob():
    ids = glob_wildcards('data/genes/BUSCO_eutheria/codon_nuc/byBUSCO/{id}.fa')
    #print(len(ids))
    #print(len(ids[0]))
    return ids[0]

rule all:
    input: 
        #expand('output/MACSE_ALFIX/BUSCO_eutheria/codon_nuc/byBUSCO/{id}.fa', id = get_BUSCO_union_ids())
        expand('output/MACSE_ALFIX/BUSCO_eutheria/codon_nuc/byBUSCO/{id}', id = get_BUSCO_union_ids_glob())

rule macse_hmmcleaner:
    input: 'data/genes/{class}/{type}/{group}/{id}.fa'
    output: 
        dir = directory('output/MACSE_ALFIX/{class}/{type}/{group}/{id}'),
        #aln = 'output/MACSE_ALFIX/{class}/{type}/{group}/{id}/{id}_final_align_NT.aln'
    log:
        out='logs/MACSE_ALFIX/{class}/{type}/{group}/{id}.log',
        err='logs/MACSE_ALFIX/{class}/{type}/{group}/{id}.err'
    shell: 
        'singularity run code/MACSE_ALFIX_v01.sif '
        ' --out_dir {output.dir} '
        ' --out_file_prefix {wildcards.id} '
        '--in_seq_file {input} '
        '--java_mem 10g 2> {log.err} > {log.out}'
