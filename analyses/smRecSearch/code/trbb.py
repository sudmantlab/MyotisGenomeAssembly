#!/usr/bin/env python3

import click
from pathlib import Path
from RecBlast.RecBlast import RecSearch

@click.command()
@click.option("-q", "--query-file", type=click.Path(exists=True))
@click.option("--query-file-type", type=str, default="fasta")
@click.option("-p", "--max-processes", type=int, default=40)
@click.option("-fp", "--forward-port")
@click.option("-rp", "--reverse-port")
@click.option("-fs", "--forward-species", type=str)
@click.option("-ft", "--forward-twobit", type=click.Path(exists=False))
@click.option("-rs", "--reverse-species", type=str)
@click.option("-rt", "--reverse-twobit", type=click.Path(exists=False))
@click.option("-ps", "--perc-score")
@click.option("-pi", "--perc-identity")
@click.option("-pq", "--perc-query-span")
@click.option("--annotation_lookup_tsv", type=str)
@click.option("--output-root", type=str, default="./output")
def __main__(query_file, forward_port, forward_species, forward_twobit,
             reverse_port, reverse_species, reverse_twobit, query_file_type="fasta",
             perc_score="0.1", perc_identity="0.5", perc_query_span="0.5", max_processes=40,
             annotation_lookup_tsv="./output/recBlastDBPrep/hg38_geneAcc_hashTable.tsv",
             output_root = "./output"):
    #perc_score = float(perc_score)
    #perc_identity = float(perc_identity)
    #perc_query_span = float(perc_query_span)

    forward_twobit = Path(forward_twobit)
    reverse_twobit = Path(reverse_twobit)

    print(forward_twobit, reverse_twobit, output_root, perc_identity, perc_score, perc_query_span, query_file, sep="\n")

if __name__ == "__main__":
    __main__()

exit()
