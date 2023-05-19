#!/usr/bin/env python3

mRNAs = [] 
exons = []
CDSs = []

with open(snakemake.input[0]) as infile:
	for line in infile:
		if line.startswith('#'):
			continue
		line = line.strip().split('\t')
		line[-1] = line[-1].split(';')
		line[-1] = (line[-1][0].replace("ID=",""), line[-1][1].replace("Parent=",""))
		if line[2] == 'mRNA':
			mRNAs.append(line[-1])
		elif line[2] == 'exon':
			exons.append(line[-1])
		elif line[2] == 'CDS':
			CDSs.append(line[-1])
		else:
			continue

CDS_mRNAs = {CDS[1] for CDS in CDSs}

missing_cds = [mRNA[0] for mRNA in mRNAs if mRNA[0] not in CDS_mRNAs]

with open(snakemake.output[0], 'w') as outfile:
	for cds in missing_cds:
		outfile.write(cds + '\n')

