from Bio import SeqIO


final_list = []
for item in SeqIO.parse(snakemake.input[0], "fasta"):
	if len(item) != 3:
		item.seq += "-"*(len(item)%3)
	final_list.append(item)


SeqIO.write(final_list, snakemake.output[0], "fasta")
