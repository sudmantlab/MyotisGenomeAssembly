#!/usr/bin/env python3.6

from pathlib import Path
import pandas as pd


def replace_column(bed_row, replace_dict, col=0):
    bed_row[col] = replace_dict[bed_row[col]]
    return bed_row


def get_replacement_dict_NCBIToUCSC(replace_table):
    rd = {rec.RefSeq_Accn: rec.GenBank_Accn.split(".")[0] for index, rec in pd.read_table(replace_table, comment="#", names=["Sequence_Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank_Accn", "Relationship", "RefSeq_Accn", "Assembly-Unit", "Sequence-Length", "UCSC_style_name"]).iterrows()}
    rd["NC_004920.1"] = "chrM"
    return rd
    
    
def get_replacement_dict_UCSCToGenBank(replace_table):
    rd = {rec.UCSC_style_name: rec.GenBank_Accn.split(".")[0]+".1" for index, rec in pd.read_table(replace_table, comment="#", names=["Sequence_Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank_Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC_style_name"]).iterrows()}
    rd["scaffold_750"] = "AAGU03094116.1"
    return rd
    
def get_replacement_dict_GenBankToUCSC(replace_table):
    rd = {rec.GenBank_Accn.split(".")[0]+".1": rec.UCSC_style_name for index, rec in pd.read_table(replace_table, comment="#", names=["Sequence_Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank_Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC_style_name"]).iterrows()}
    #rd["chrM"] = rd["MT"]
    return rd

def __main__(infile, outfile, replace_table, frm="UCSC_style_name"):
    inpath = Path(infile)
    outpath = Path(outfile)
    replace_dict = get_replacement_dict_UCSCToGenBank(replace_table) if frm == "UCSC_style_name" else get_replacement_dict_GenBankToUCSC(replace_table)
    with inpath.open() as bed, outpath.open("w") as outf:
        l = ("\t".join(replace_column(line.strip().split("\t"), replace_dict, col=0)) for line in bed if "chrM" not in line)
        outf.write("\n".join(l))


if __name__ == "__main__":
    __main__(infile=snakemake.input["bed"], outfile=snakemake.output[0], replace_table=snakemake.input["replace_table"], frm=snakemake.params["frm"])
