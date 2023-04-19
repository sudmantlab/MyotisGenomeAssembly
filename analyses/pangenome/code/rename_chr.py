#!/usr/env python3

from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd

#genome_dict = SeqIO.index(snakemake.input["genome"], "fasta")
genome_list = [rec for rec in SeqIO.parse(snakemake.input["genome"], "fasta")]

flip_table = pd.read_csv(snakemake.input["flip_table"], comment="#")

print(flip_table)

def rename_chr(d):
    # d should be a dictionary that packages the scaffold, and the flip_table
    # pangenome naming format: genome#hap#scaffold#og_scaffold
    # the reference genome is tname
    # the genome of interest is qname
    
    #print(d['chr'].id + " start")
    
    g = snakemake.wildcards.genome
    hap = "0"
    #print(d['chr'].id + " tname")
    og_name = d['chr'].id
    ref_name =  d['tbl'][d['tbl'].qname == og_name]['tname'].values[0]
    new_id = "#".join((g, hap, ref_name, og_name))
    
    #print(d['chr'].id + " flip")
    do_flip = d['tbl'][d['tbl'].qname == og_name]['flip'].values[0]
    
    new_chr = d['chr'].reverse_complement() if do_flip else d['chr']
    new_chr.description = 'flipped' if do_flip else ''
    new_chr.id = new_id
    
    return new_chr.format("fasta")


if __name__ == '__main__':
    dict_args_list = []
    for rec in genome_list:
        if (flip_table['qname'].eq(rec.id)).any():
            dict_args_list.append(dict(tbl = flip_table, chr = rec))
    
    with Pool(processes=snakemake.threads) as pool:
        new_genome = pool.map(rename_chr, dict_args_list)
    with open(snakemake.output[0], 'w') as out_genome:
       out_genome.writelines(new_genome)
    
