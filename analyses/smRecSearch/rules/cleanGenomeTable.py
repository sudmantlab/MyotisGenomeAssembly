#!/usr/bin/env python3.6

import pandas as db
from numpy import nan

def ucsc_name_gen(species_name):
    parts = [s[0:3] for s in species_name.lower().split("_")] + ["1"]
    parts[1] = parts[1].capitalize()
    return "".join(parts)

def __main__(hal_db, ncbi_db, outfile):
  db1 = pd.read_table(hal_db, names=['species','AssemblyAccession'])
  db1_discovar = db1[db1.AssemblyAccession == "DISCOVAR"].assign(Synonym__Genbank=nan, UCSCName=nan)
  db2 = pd.read_csv(ncbi_db, usecols=['AssemblyAccession', 'Synonym__Genbank', 'SpeciesName', "UCSCName"]).rename(columns={"SpeciesName": 'species'}).assign(species = lambda df: [spc.replace(" ", "_") for spc in df.species])
  dbA = db1.merge(db2, how="inner", left_on="AssemblyAccession", right_on="Synonym__Genbank", suffixes=("_hal", "_ncbi"))[['species_hal', 'AssemblyAccession_hal', 'UCSCName', 'species_ncbi', 'Synonym__Genbank']].rename(columns={'AssemblyAccession_hal':'AssemblyAccession'})
  dbB = db1.merge(db2, how="inner", left_on="AssemblyAccession", right_on="AssemblyAccession", suffixes=("_hal", "_ncbi"))
  db3 = pd.concat([pd.concat([dbA, dbB]).drop_duplicates(keep=False)[['species_hal', 'AssemblyAccession', 'UCSCName']].rename(columns={"species_hal":"species"}), db1_discovar]).assign(UCSCName = lambda df: [ucsc_name_gen(row[0]) if pd.isna(row[1]) else row[1] for row in zip(df.species, df.UCSCName)])
  db3.to_csv(outfile)
  

if __name__ == "__main__":
    __main__(snakemake.input["hal_db"], snakemake.input("ncbi_db"), snakemake.output[0])
