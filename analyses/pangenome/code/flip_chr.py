#!/usr/env python3

from multiprocessing import Pool
from Bio import SeqIO
import pandas as pd

genome_dict = SeqIO.index(snakemake.input["genome"], "fasta")

flip_table = pd.read_csv(snakemake.input["flip_table"])

print(flip_table)

def flip_chr(d):
    new_chr = d['chr'].reverse_complement() if d['do_flip'] else d['chr']
    new_chr.id = d['chr'].id
    new_chr.description = 'flipped' if d['do_flip'] else ''
    return new_chr.format("fasta")

#with open(snakemake.output[0], 'w') as out_genome:
#    for row in zip(flip_table['qname'], flip_table['flip']): 
#        print(row)
#        print(row[0])
#        print(row[1])
#        new_chr = genome_dict[row[0]].reverse_complement() if row[1] else genome_dict[row[0]]
#        out_genome.write(new_chr.format("fasta"))

if __name__ == '__main__':
    with Pool(processes=snakemake.threads) as pool:
        new_genome = pool.map(flip_chr, (dict(chr=genome_dict[row[0]], do_flip=row[1]) for row in zip(flip_table['qname'], flip_table['flip'])))
    with open(snakemake.output[0], 'w') as out_genome:
       out_genome.writelines(new_genome)
    
