#!/usr/bin/env python3

import click
import sys
from pathlib import Path
from RecBlast.RecBlast import RecSearch
import RecBlast.WarningsExceptions as RBWE


def deduce_searchtype(query_type, db_type, search_algorithm):
    # a bit of cleaning
    query_type = query_type.lower()
    db_type = db_type.lower()
    search_algorithm = search_algorithm.lower()

    if "blast" in search_algorithm:
        if query_type == "dna":
            if db_type == "prot":
                return "blastx"
            elif db_type == "dna":
                return "blastn"
            else:
                raise Exception("Unknown search database type! Allowed options are 'dna' or 'prot'")
        elif query_type == "prot":
            if db_type == "prot":
                return "blastp"
            elif db_type == "dna":
                return "tblastn"
            else:
                raise Exception("Unknown search database type! Allowed options are 'dna' or 'prot'")
        else:
                raise Exception("Unknown search sequence type! Allowed options are 'dna' or 'prot'")
    if "blat" in search_algorithm:
        if query_type == "dna":
            if db_type == "prot":
                return "blatx"
            elif db_type == "dna":
                return "blat"
            else:
                raise Exception("Unknown search database type! Allowed options are 'dna' or 'prot'")
        elif query_type == "prot":
            if db_type == "prot":
                return "blatp"
            elif db_type == "dna":
                return "tblat"
            else:
                raise Exception("Unknown search database type! Allowed options are 'dna' or 'prot'")
        else:
                raise Exception("Unknown search sequence type! Allowed options are 'dna' or 'prot'")
    else:
        raise RBWE.SearchEngineNotImplementedError("This search engine hasn't been implemented yet! Only BLAT and BLAST have been implemented!")


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
@click.option("-ps", "--perc-score", type=str, default= "0.1")
@click.option("-pi", "--perc-identity", type=str, default = "0.5")
@click.option("-pq", "--perc-query-span", type=str, default = "0.5")
@click.option("--query_type", type=str, default = "prot")
@click.option("--reverse_type", type=str, default = "dna")
@click.option("--forward_algo", type=str, default = "blat")
@click.option("--reverse_algo", type=str, default = "blat")
@click.option("--reverse_db_type", type=str, default = "dna")
@click.option("--forward_db_type", type=str, default = "dna")
@click.option("--annotation_lookup_tsv", type=str, default = "")
@click.option("--output-root", type=str, default="./output")
@click.option('-v', '--verbose', count=True)
def __main__(query_file, forward_port, forward_species, forward_twobit,
             reverse_port, reverse_species, reverse_twobit, 
             query_type, forward_db_type, forward_algo, 
             reverse_type, reverse_db_type, reverse_algo,
             perc_score, perc_identity, perc_query_span, query_file_type, max_processes,
             annotation_lookup_tsv, output_root, verbose):
    perc_score = float(perc_score)
    perc_identity = float(perc_identity)
    perc_query_span = float(perc_query_span)

    forward_twobit = Path(forward_twobit)
    reverse_twobit = Path(reverse_twobit)

    print(forward_twobit, reverse_twobit, output_root, perc_identity, perc_score, perc_query_span, query_file, sep="\n", file=sys.stderr)
    output_location = Path(output_root, forward_twobit.stem)
    print(output_location, file=sys.stderr)
    
    f_search_type = deduce_searchtype(query_type, forward_db_type, forward_algo)
    r_search_type = deduce_searchtype(reverse_type, reverse_db_type, reverse_algo)

    recblast = RecSearch(target_species=forward_species, query_species=reverse_species,
                         forward_search_type=f_search_type, reverse_search_type=r_search_type,
                         sequence_source="twobit", verbose=verbose)
    recblast.max_processes = max_processes
    recblast.set_queries(query_file,
                         infile_type=query_file_type)
    recblast.forward_search_settings['database_port'] = {forward_species: forward_port}
    recblast.forward_search_settings['database'] = {forward_species: str(forward_twobit.name)}
    recblast.forward_search_settings['database_path'] = str(forward_twobit.parent)
    recblast.forward_search_criteria = dict(perc_score=perc_score,
                                            perc_ident=perc_identity,
                                            perc_query_span=perc_query_span)
    recblast.sequence_source_settings['database'] = {forward_species: str(forward_twobit.name)}
    recblast.sequence_source_settings['database_path'] = str(forward_twobit.parent)
    recblast.memory_saver_level = 1
    recblast.reverse_search_settings['database'] = {reverse_species: str(reverse_twobit.name)}
    recblast.reverse_search_settings['database_path'] = str(reverse_twobit.parent)
    recblast.reverse_search_settings['database_port'] = {reverse_species: reverse_port}
    if annotation_lookup_tsv:
        recblast.set_translation_annotation_parameters(method="table", key_value_order=False,
                                                    tsv_location=annotation_lookup_tsv)
    else:
        recblast.set_translation_annotation_parameters(method=False)
    recblast(run_name="{0}-pcScore{1}_pcIdent{2}_pcQuerySpan{3}_reverse-{4}".format(Path(query_file).stem, 
                                                                                    perc_score, 
                                                                                    perc_identity, 
                                                                                    perc_query_span, 
                                                                                    reverse_twobit.stem),
             output_type="bed-complete",
             output_location=output_location)
    

if __name__ == "__main__":
    __main__()

exit()
