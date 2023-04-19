#!/usr/bin/env python3

import click
from pathlib import Path


@click.command()
@click.argument('files', nargs=-1)
def __main__(files):
    for rbb_file in files:
        inpath = Path(rbb_file)
        outpath = inpath.parent.joinpath("..", "RBB", inpath.name + ".rbb")
        droppath = inpath.parent.joinpath("..", "RBB", "dropped", inpath.name + ".dropped")
        try:
            outpath.parent.mkdir(parents=True)
        except FileExistsError:
            pass
        try:
            droppath.parent.mkdir(parents=True)
        except FileExistsError:
            pass
        with inpath.open() as bed, outpath.open("w") as outf, droppath.open("w") as dropped:
            uniq_rec = {}
            for line in bed:
                try:
                    record = line.strip().split("\t")
                    bedline = record[0:12]
                    chr = record[0]
                    s = record[1]
                    e = record[2]
                    name = record[3]
                    query = record[3].split("_")[0]
                    score = record[4]
                    strand = record[5]
                    thickStart = record[6]
                    thickEnd = record[7]
                    itemRbg = record[8]
                    blockCount = record[9]
                    blockSizes = record[10]
                    blockStarts = record[11]

                    if "annotations" in record[-1]:
                        #print(record[-1])
                        anno = record[-1].split("annotations=")[1].strip()
                        if query in anno:
                            anno_set = [chr for chr in anno.split(",") if query in chr][0]
                            #print(anno_set)
                            if ":" in anno_set:
                                anno_chr = anno_set.split(":")[0]
                                #print(anno_chr)
                                anno_s = anno_set.split(":")[1].split("-")[0]
                                #print(anno_s)
                                anno_e = anno_set.split(":")[1].split("-")[1]
                            else:
                                anno_chr = anno_set
                                anno_s = ""
                                anno_e = ""
                            #print(anno_e)
                            bedline += [anno_chr, anno_s, anno_e]
                        else:
                            #print(type(anno))
                            dropped.write("\t".join(bedline) + "\n")
                            continue
                    else:
                        #print(type(anno))
                        dropped.write("\t".join(bedline) + "\n")
                        continue
                    if "query_coverage" in record[-1]:
                        import re
                        p = re.compile("(query_coverage=(\(\d+, \d+\),?)+)")
                        qc_ranges = p.match(record[-1]).group()
                        bedline.append(qc_ranges)
                    outf.write("\t".join(bedline) + "\n")
                except IndexError as err:
                    print(err)
                    dropped.write(line)

if __name__ == "__main__":
    try:
        __main__(*snakemake.input)
    except NameError:
        __main__()

