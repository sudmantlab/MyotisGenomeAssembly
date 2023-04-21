from Bio import SeqIO

f_infasta = snakemake.input['fasta']
f_tt = snakemake.input['tt']
f_outfile = snakemake.output[0]
genome = snakemake.wildcards['genome']

with open(f_tt) as infile:
	tt = {v:k for k, v in (i.strip().split('\t') for i in infile)}

records = list(SeqIO.parse(f_infasta, 'fasta'))

with open(f_outfile, 'w') as outfile:
	for record in records:
		record.description = record.id
		record.id = genome + "_" + tt[record.id]
	SeqIO.write(records, outfile, 'fasta-2line')


