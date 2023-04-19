#outfile = 'output/contigs_to_drop/repeatMasker_denovo/mMyoVel1.repeat_over_90pc'
#genomefile = 'data/genomes/mMyoVel1.genome'
#rm_gff = '../data/repeatMasker/denovo/mMyoVel1.denovo.out.reformated.gff'

genomefile = snakemake.input.genome
rm_gff = snakemake.input.gff
outfile = snakemake.output[0]

list_contig = []

with open(genomefile) as genome, open(rm_gff) as gff:
    g = [i.strip().split() for i in genome]
    f = [i.strip().split() for i in gff]

for line in f:
    for contig, l in g:
        if line[0] != contig:
            continue
        len = int(line[4]) - int(line[3])
        if len >= 0.9*int(l):
            list_contig.append(contig)

with open(outfile,'w') as outf:
    outf.write('\n'.join(list_contig))


