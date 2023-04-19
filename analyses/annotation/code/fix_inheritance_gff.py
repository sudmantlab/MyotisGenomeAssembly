#!/usr/bin/env python3

from collections import Counter
from copy import copy

id_count = Counter()
target_count = Counter() #just in case
mp_count = Counter() #also just in case
with open(snakemake.input[0]) as infile, open(snakemake.output[0], "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            continue
        line = line.strip()
        cols = line.split('\t')
        attr = cols[-1].split(';')
        if cols[2] == "mRNA":
            id = attr[0].replace("ID=","")
            mp_count[id] += 1
            #if mp_count >1:
            #    id += "_" + str(mp_count[id])
            target = attr[-1].replace(' ', '_').replace('Target=','')
            target_count[target] += 1
            if target_count[target] >1:
                target += '_' + str(target_count[target])
            geneattr = copy(attr)
            geneattr[0] = "ID=" + target
            genecols = copy(cols)
            genecols[2] = "gene"
            genecols[-1] = ";".join(geneattr)
            geneline = "\t".join(genecols) + "\n"
            attr[0] = 'ID=' + id + ";Parent=" + target
            outfile.write(geneline)
        else:
            id = attr[0].replace("Parent=", "")
            id_count[id] += 1
            attr[0] = attr[0] + ";ID=" +  id + "." + str(id_count[id])
            exoncols = copy(cols)
            exonattr = copy(attr)
            exoncols[2] = 'exon'
            exoncols[7] = '.'
            exonattr[0] += '.exon'
            exoncols[-1] = ';'.join(exonattr)
            exonline = '\t'.join(exoncols) + '\n'
            outfile.write(exonline)
        cols[-1] = ";".join(attr)
        newline = "\t".join(cols) + "\n"
        outfile.write(newline)
