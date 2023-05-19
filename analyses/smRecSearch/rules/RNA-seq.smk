import pandas as pd
import numpy

########################################################################################################################
######## Config
########################################################################################################################





def genome_to_species(genome, genome_df = pt):
    selection = genome_df['Genome'] == genome
    return genome_df.loc[selection]["Species"][0]


def get_SRA_list(sra_file):
    """Reads a .csv SRA file and returns a dataframe indexed by SRA Run accessions"""
    df = pd.read_csv(sra_file, usecols=['Run', 'Organism'], comment='#').set_index("Run", drop=False)
    return df


def get_SRA_acc(species, dir_sra, suffix="_SraRunTable.tsv"):
    """Returns a list of SRA Run Accession numbers from a SraRunTable"""
    try:
        species = species.replace(" ", "_")
        sra_table = species + suffix
        sra_file = dir_sra.rstrip("/") + "/" + sra_table
        df = get_SRA_list(sra_file)
        return list(df.index)
    except FileNotFoundError:
        return ""
    except AttributeError:
        return ""

def get_SRA_acc2(species, sra_file= config["sra_file"]):
    """Returns a list of SRA Run Accession numbers from a master SraRunTable file"""
    df = get_SRA_list(sra_file)
    selection = df["Organism"] == species
    sra_acc_list = df[selection]
    return list(sra_acc_list.index)


def generate_genome_SRA_list(genome_df, sra_file= config["sra_file"]):
    """Generates a set of species-SRA strings from a list of species"""
    genome_sra="{genome}_{sra}"
    final_list = []
    for species, genome in zip(genome_df.index, genome_df):
        # sra_list = get_SRA_acc(species, dir_sra, suffix)
        sra_list = get_SRA_acc2(species, sra_file)
        if sra_list:
            final_list.extend(map(lambda sra: (genome, sra), sra_list))
        else:
            continue
    if final_list == []:
        raise Exception("generate_genome_SRA_list error: There's no SRAs for any of the provided species!")
    return [genome_sra.format(genome=i[0], sra=i[1]) for i in final_list]

def sra_from_wildcards(wildcards, sra_file= config["sra_file"], genome_df = pt):
    selection = genome_df['Genome'] == wildcards.genome
    return get_SRA_acc2(genome_df.loc[selection].index[0], sra_file)

def species_hasSRA(genome, sra_file=config["sra_file"], genome_df = pt):
    species = genome_to_species(genome)
    df = pd.read_csv(sra_file, usecols=['Organism'], comment='#')
    return (species == df.Organism).any()


def stmerge_inputs(wildcards, sra_file= config["sra_file"], genome_df = pt):
    sras = sra_from_wildcards(wildcards, genome_df = pt, sra_file= config["sra_file"])
    final_list = ["{path}/{genome}_{sra}.denovo.gff3".format(path = dir_gff3,
                                                genome=wildcards.genome,
                                                sra=s)
                for s in sras
           ]
    if final_list == []:
        raise Exception("STMerge-Denovo Error: There's no SRAs for the selected species!")
    else:
        return final_list

def stmerge_final_inputs(wildcards, sra_file = config["sra_file"], genome_df = pt):
    sras = sra_from_wildcards(wildcards, genome_df = pt, sra_file= config["sra_file"])
    final_list = ["{path}/{genome}_{sra}-final.gff3".format(path = dir_stringtie,
                                                genome=wildcards.genome,
                                                sra=s)
                for s in sras
           ]
    if final_list == []:
        raise Exception("STMerge-Final Error: There's no SRAs for the selected species!")
    else:
        return final_list

loxAfr_genomes = pt[pt.Genome.isin(["loxAfr3","loxAfr4"])].Genome

########################################################################################################################
######## Rules
########################################################################################################################

rule all_alignments:
    input:
        expand("data/BED/{genome}-finalGuide.bed", genome=genomes)

rule loxAfr_alignments:
    """Useful for troubleshooting"""
    input:
        expand("data/BED/{genome}-finalGuide.bed", genome=loxAfr_genomes)
#    output:
#       "./loxAfr_alignments.list"
#    shell:
#        "echo {input} >> {output}"


rule hisat2_index:
    input:
        "{path}/{genome}{extension}".format(genome='{genome}', path=dir_genome, extension = ".fa")
    output:
        temp("{path}/{genome}{extension}".format(genome='{genome}', path=dir_idx, extension = ".{dataset,\d+}.ht2"))
    params:
        real_output="{path}/{genome}{extension}".format(genome='{genome}', path=dir_idx, extension = "")
    threads: 20
    conda: "../envs/conda_tuxedo.yaml"
    log: "{path}/{genome}{extension}".format(genome='{genome}', path=log_idx, extension = "{dataset,\d+}.log")
    shell:
        "hisat2-build -p {threads} {input} {params.real_output}"

rule hisat2:
    input:
        idx_fake = "{path}/{genome}{extension}".format(genome='{genome}', path=dir_idx, extension = ".1.ht2")
    output:
        temp("{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=dir_bam, extension = ".bam"))
    log:
        "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=log_bam, extension = ".bam.log")
    params:
        summary = "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=log_bam, extension = ".bam.summary"),
        idx = "{path}/{genome}".format(genome='{genome}', path=dir_idx)
    threads: 20
    conda: "../envs/conda_tuxedo.yaml"
    shell:
        "hisat2 -p {threads} "
        "-x {params.idx} "
        "--dta "
        "--summary-file {params.summary} "
        "--met-stderr "
        "--sra-acc {wildcards.sra} "
        "2> {log} | "
        "samtools view -b -u | "
        "samtools sort - "
        "--threads {threads} "
        "-o {output}"

rule samtools_index:
    input:
        bam = "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=dir_bam, extension = ".bam")
    output:
        bai = temp("{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=dir_bam, extension = ".bai"))
    log: "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=log_bam, extension = ".bai.log")
    threads: 10
    conda: "../envs/conda_tuxedo.yaml"
    shell:
        "samtools index {input.bam}"

rule stringtie_denovo:
    input:
        "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=dir_bam, extension = ".bam")
    output:
        temp("{path}/{genome}_{sra}.denovo{extension}".format(genome='{genome}', sra='{sra}', path=dir_gff3, extension = ".gff3"))
    log: "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=log_gff3, extension = ".gff3.log")
    params:
    threads: 10
    conda: "../envs/conda_tuxedo.yaml"
    shell:
        "stringtie {input} "
        "-p {threads} "
        "-l {wildcards.genome}_{wildcards.sra} "
        "-v "
        "-o {output} "
        "2> {log}"

rule stringtie_merge:
    input: stmerge_inputs
    output: "{path}/{genome}{extension}".format(genome='{genome}', path=dir_stringtieMerge, extension = "-guide.gff3")
    log: "{path}/{genome}{extension}".format(genome='{genome}', path=log_stringtie, extension = "-guide.gff3.log")
    threads: 10
    conda: "../envs/conda_tuxedo.yaml"
    shell:
        "stringtie --merge "
        "{input} "
        "-l {wildcards.genome} "
        "-v "
        "-o {output} "
        "2> {log}"

rule stringtie_final:
    input:
        bam = "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=dir_bam, extension = ".bam"),
        guide = "{path}/{genome}{extension}".format(genome='{genome}', path=dir_stringtieMerge, extension = "-guide.gff3")
    output:
        protected("{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=dir_stringtie, extension = "-final.gff3"))
    log: "{path}/{genome}_{sra}{extension}".format(genome='{genome}', sra='{sra}', path=log_stringtie, extension = "-final.gff3.log")
    conda: "../envs/conda_tuxedo.yaml"
    threads: 10
    shell:
        "stringtie {input.bam} "
        "-p {threads} "
        "-G {input.guide} "
        "-l {wildcards.genome} "
        "-vBAC "
        "-o {output} "
        "2> {log}"

rule stringtie_merge_final:
    input: stmerge_final_inputs
    output: "{path}/{genome}{extension}".format(genome='{genome}', path=dir_stringtieMergeFinal, extension = "-finalGuide.gff3")
    log: "{path}/{genome}{extension}".format(genome='{genome}', path=log_stringtieMergeFinal, extension = "-finalGuide.gff3.log")
    threads: 10
    conda: "../envs/conda_tuxedo.yaml"
    shell:
        "stringtie --merge "
        "{input} "
        "-l {wildcards.genome} "
        "-v "
        "-o {output} "
        "2> {log}"


rule GFF3ToBed:
    input:
        "{path}/{genome}{extension}".format(genome='{genome}', path=dir_stringtieMergeFinal, extension = "-finalGuide.gff3")
    output:
        "output/BED/{genome}{extension}".format(file="{file}", genome="{genome}", extension = "-finalGuide.bed")
    conda: "../envs/conda_tuxedo.yaml"
    shell:
        "gff2bed < {input} > {output}"

rule linkBed:
    input: "output/BED/{file}.bed"
    output: "data/BED/{file}.bed"
    shell: "ln -sf ../../{input} {output}"
