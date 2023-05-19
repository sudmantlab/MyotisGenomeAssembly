import vcf
import sys

infile_name = sys.argv[1]
outfile_name = sys.argv[2]

print("<", infile_name, file=sys.stderr)
print(">", outfile_name, file=sys.stderr)

def format_het(record):
    line = [
            record.CHROM,
            record.POS-1,
            record.POS,
            record.heterozygosity
           ]
    return "\t".join((str(i) for i in line))+"\n"

with open(infile_name) as infile, open(outfile_name,"w") as outfile:
    vcf_reader = vcf.Reader(infile)
    for rec in vcf_reader:
        if rec.heterozygosity != 0.0:
            outfile.write(format_het(rec))


