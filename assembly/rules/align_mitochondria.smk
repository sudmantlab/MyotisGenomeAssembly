def get_all_mitoHiFi_genomes(wildcards):
    samples = pd.read_table("pep_mitochondrion.tsv", index_col="species")
    return samples.loc[wildcards.species].
