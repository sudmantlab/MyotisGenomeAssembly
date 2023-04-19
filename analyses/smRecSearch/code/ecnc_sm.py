#!/usr/bin/env python3.6

import re
from pathlib import Path
from RecBlast import merge_ranges


def ecnc(ranges):
    coverage = merge_ranges(ranges)
    sum_coverage = sum([i[1] - i[0] for i in coverage])
    sum_nuc = sum([r[1] - r[0] for r in ranges])
    return round(sum_nuc/sum_coverage, 2)

def __main__(infile, outfile, dropfile):
    p = re.compile("(?<=query_coverage=)(\(\d+, \d+\),?)+")
    newp = False
    inpath = Path(infile)
    outpath = Path(outfile)
    droppath = Path(dropfile)
    try:
        outpath.parent.mkdir(parents=True)
    except FileExistsError:
        pass
    try:
        droppath.parent.mkdir(parents=True)
    except FileExistsError:
        pass
    with inpath.open() as bed, outpath.open("w") as outf, droppath.open("w") as dropped:
        hit_dict = {}
        for line in bed:
            try:
                record = line.strip().split("\t")
                name = record[3]
                query = name.split("_")[0]
                qcov = p.search(record[-1])
                if qcov is None:
                    p = re.compile("(?<=query_coverage=)(\(\d+,\d+\),?)+")
                    qcov = p.search(record[-1])
                    newp = True
                
                if newp: 
                    qc_ranges = ((int(i.split(',')[0].lstrip("(")), int(i.split(',')[1].rstrip(")"))) for i in
                                 qcov.group().strip("()").split("),("))
                else:
                    qc_ranges = ((int(i.split(', ')[0].lstrip("(")), int(i.split(', ')[1].rstrip(")"))) for i in
                                 qcov.group().strip("()").split("),("))
                qc_ranges_sorted = [sorted(i) for i in qc_ranges]
                if query in hit_dict:
                    hit_dict[query] += qc_ranges_sorted
                else:
                    hit_dict[query] = qc_ranges_sorted
            
            except AttributeError as err:
                print(err)
                print(line)
                dropped.write(line)
            except Exception as err:
                print(err)
                print(line)
                dropped.write(line)
        hit_ecnc = ("\t".join((q, str(ecnc(r)))) for q, r in hit_dict.items())
        outf.write("\n".join(hit_ecnc))


if __name__ == "__main__":
    __main__(snakemake.input[0], snakemake.output["ECNC_file"], snakemake.output["dropped_lines"])
