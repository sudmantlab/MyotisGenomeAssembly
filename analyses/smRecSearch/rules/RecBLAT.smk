import pandas as pd
from snakemake.exceptions import WorkflowError

def get_port_for_genome(df, genome, translated = True):
    try:
        if translated:
            return int(df.loc[df["Genome"]==genome].Translated_BLAT_Port)
        else:
            return int(df.loc[df["Genome"]==genome].Untranslated_BLAT_Port)
    except TypeError:
        raise Exception("Port number for genome {} not found; check your spelling or add the genome to the portTable file!".format(genome))

def get_species_port_translated(df):
    return ["{{path}}/translated-{0}-{1}".format(i, j) for i, j in zip(df.Translated_BLAT_Port, df.Genome)]

def get_species_port_untranslated(df):
    return ["{{path}}/untranslated-{0}-{1}".format(i, j) for i, j in zip(df.Untranslated_BLAT_Port, df.Genome)]

def get_species_port_closed(df):
    l = []
    l += ["{{path}}/translated-{0}-{1}-closed".format(j, i) for i, j in zip(df.Genome, df.Translated_BLAT_Port)]
    l += ["{{path}}/untranslated-{0}-{1}-closed".format(j, i) for i, j in zip(df.Genome, df.Untranslated_BLAT_Port)]
    return l

def get_port(f):
    with open(f) as infile:
        return int(infile.read())

#def filtered_beds(wildcards):
#    checkpoint_output = checkpoints.filterRBHBbyBED.get(**wildcards).output[0]
#    return expand("output/RBB/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected-{{type}}.bed", 
#                  genome = wildcards.genome, query = wildcards.query, pc_score = wildcards.pc_score, pc_ident=wildcards.pc_ident, pc_qspan = wildcards.pc_qspan, rgenome = wildcards.rgenome, 
#                  type = glob_wildcards(checkpoint_output).type)

def filtered_beds(wildcards):
    t = glob_wildcards("data/BED/{genome}-{{type}}.bed".format(genome=wildcards.genome)).type
    has_SRA = species_hasSRA(wildcards.genome)
    if t:
        final_list = expand("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected-{type}.bed",
                      genome = wildcards.genome, query = wildcards.query, pc_score = wildcards.pc_score, pc_ident=wildcards.pc_ident, pc_qspan = wildcards.pc_qspan, rgenome = wildcards.rgenome,
                      type = t)
    else:
        final_list = expand("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb",
                      genome = wildcards.genome, query = wildcards.query, pc_score = wildcards.pc_score, pc_ident=wildcards.pc_ident, pc_qspan = wildcards.pc_qspan, rgenome = wildcards.rgenome)
    if has_SRA:
        final_guide = expand("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected-finalGuide.bed",
                      genome = wildcards.genome, query = wildcards.query, pc_score = wildcards.pc_score, pc_ident=wildcards.pc_ident, pc_qspan = wildcards.pc_qspan, rgenome = wildcards.rgenome
                      )
    else:
        final_guide = []
    return final_list + final_guide


def translated_BLAT_db_type(query_type, db_type):
    """ Only returns true if the query is Protein and DB is DNA. Note that in the case of query=DNA and db=prot, gfClient is the one that has to implement dnax"""
    # a bit of cleaning
    query_type = query_type.lower()
    db_type = db_type.lower()
    if query_type == "dna":
        if db_type == "prot":
            return False
        elif db_type == "dna":
            return False
        else:
            raise Exception("Unknown search database type! Allowed options are 'dna' or 'prot'")
    elif query_type == "prot":
        if db_type == "prot":
            return False
        elif db_type == "dna":
            return True
        else:
            raise Exception("Unknown search database type! Allowed options are 'dna' or 'prot'")
    else:
        raise Exception("Unknown search sequence type! Allowed options are 'dna' or 'prot'")
# Parameters

# Output Directories


# Log Directories

rule RBHB_loxAfr3:
    input:
        expand("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_evidenced.bed",
               genome="loxAfr3", query="AvA", pc_score = 0.1, pc_ident = 0.8, pc_qspan = 0.5, rgenome = "hg38_maskRep_noVarChr_fragWithGenes")

rule start_translated_gfServer:
    input:
        twoBitFile=ancient("{path}/{{genome}}{extension}".format(path=dir_2bit, extension = ".2bit"))
    output:
        "{path}/translated-{{genome}}".format(path=dir_flags)
    params:
        twobitdir = os.getcwd() + "/" + dir_2bit
    log:
        "{path}/translated-{{genome}}.log".format(path=log_flags)
    group: "recBlat"
    shell:
        "./code/start_Translated_gfServer.sh {wildcards.genome} {output} {params.twobitdir} > {log}"

rule start_untranslated_gfServer:
    input:
        twoBitFile=ancient("{path}/{{genome}}{extension}".format(path=dir_2bit, extension = ".2bit"))
    output:
        "{path}/untranslated-{{genome}}".format(path=dir_flags)
    log:
        "{path}/untranslated-{{genome}}.log".format(path=log_flags)
    params:
        twobitdir = os.getcwd() + "/" + dir_2bit
    group: "recBlat"
    shell:
        "./code/start_Untranslated_gfServer.sh {wildcards.genome} {output} {params.twobitdir} > {log}"

def get_portfile_forward(wildcards):
    trans = translated_BLAT_db_type(config["query_type"], config["forward_db_type"])
    if trans:
        return "{path}/translated-{genome}".format(path=dir_flags, genome = wildcards.genome)
    else:
        return "{path}/untranslated-{genome}".format(path=dir_flags, genome = wildcards.genome)

def get_portfile_reverse(wildcards):
    trans = translated_BLAT_db_type(config["forward_db_type"], config["reverse_db_type"])
    if trans:
        return "{path}/translated-{genome}".format(path=dir_flags, genome = wildcards.rgenome)
    else:
        return "{path}/untranslated-{genome}".format(path=dir_flags, genome = wildcards.rgenome)


rule RecBlat:
    input:
        qfile = ancient("data/input/{query}.fa"),
        forwardTwoBitFile = "data/2bit/{genome}.2bit",
        reverseTwoBitFile = "data/2bit/{rgenome}.2bit",
        tportfile = get_portfile_forward,
        utportfile = get_portfile_reverse
    output:
        protected("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/output/{genome}_RecBlastOutput.bed")
    params:
        tport = lambda wildcards, df=pt: get_port_for_genome(df, 
                                                             wildcards.genome, 
                                                             translated=translated_BLAT_db_type(config["query_type"], 
                                                                                                config["forward_db_type"])),
        utport = lambda wildcards, df=pt: get_port_for_genome(df, 
                                                              wildcards.rgenome, 
                                                              translated=translated_BLAT_db_type(config["forward_db_type"], 
                                                                                                 config["reverse_db_type"])),
        query_type = config["query_type"],
        reverse_type = config["reverse_type"],
        forward_algo = config["forward_algo"],
        reverse_algo = config["reverse_algo"],
        reverse_db_type = config["reverse_db_type"],
        forward_db_type = config["forward_db_type"],
        annoTable = "".join(("--annotation_lookup_tsv ", config["anno_table"])) if config["anno_table"] else "",
        verbose = config["verbose"] if config["verbose"] else "v"
    threads: 40
    log: "logs/RBB/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/output/{genome}_RecBlastOutput.log"
    group: "recBlat"
    conda: "../envs/conda_recBlat.yaml"
    shell:
        "(./code/rbb.py "
        "--query-file {input.qfile} "
        "--forward-port {params.tport} "
        "--reverse-port {params.utport} "
        "--forward-species {wildcards.genome} "
        "--forward-twobit {input.forwardTwoBitFile} "
        "--reverse-species 'Homo sapiens' "
        "--reverse-twobit {input.reverseTwoBitFile} "
        "{params.annoTable} "
        "--perc-score {wildcards.pc_score} "
        "--query_type {params.query_type} "
        "--reverse_type {params.reverse_type} "
        "--forward_algo {params.forward_algo} "
        "--reverse_algo {params.reverse_algo} "
        "--forward_db_type {params.forward_db_type} "
        "--reverse_db_type {params.reverse_db_type} "
        "--perc-identity {wildcards.pc_ident} "
        "--perc-query-span {wildcards.pc_qspan} "
        "--max-processes {threads}; "
        "sed -i.bak 's/, /,/g' {output}"
        ") 2> {log}"

rule translateBedNames:
    input: 
        bed = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/output/{genome}_RecBlastOutput.bed",
        transtable = "data/input/uniprot-myoLuc-transtable.tsv"
    output: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/output_clean/{genome}_RecBlastOutput.bed"
    run:
        import sys
        with open(input.transtable) as tt:
            tdict = {l[0]: l[1] for l in [line.strip().split("\t") for line in tt]}

        with open(input.bed) as inbed, open(output[0], "w") as outbed:
            for line in inbed:
                line = line.strip().split("\t")
                id_full = line[3].split("_")
                if id_full[1] in tdict.keys():
                    line[3] = tdict[id_full[1]] + "_" + id_full[-1]
                else:
                    error = "Item " + line[3] + " is missing in translation table!"
                    print(error, file=sys.stderr)
                line = "\t".join(line) + "\n"
                outbed.write(line)

rule reciprocalBestHits:
    input:
        ancient("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/output_clean/{genome}_RecBlastOutput.bed")
#        "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/output/{genome}_RecBlastOutput.bed"
    output:
        "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb"
    conda: "../envs/conda_recBlat.yaml"
    shell:
        "./code/rbhb_from_bed.py {input} && sed -i.bak 's/, /,/g' {output}"

rule filterRBHBbyBED:
    input:
        a=ancient("output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb"),
        b=ancient("data/BED/{genome}-{type}.bed")
    output:
        "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected-{type}.bed"
    log:
        "logs/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected_{type}.log"
    conda: "../envs/conda_recBlat.yaml"
    shell:
        "bedtools intersect -u -a {input.a} -b {input.b} | bedtools sort | sed 's/, /,/g' > {output} 2>{log}"

rule evidencedRBHB:
    input:
        ancient(filtered_beds)
    output:
        "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_evidenced.bed"
    log:
        "logs/RBB/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_evidenced.log"
    conda: "../envs/conda_recBlat.yaml"
    shell:
        "cat {input} | sort | uniq | bedtools sort | sed 's/, /,/g' > {output}"
        #"echo {input} > {output}"

rule ECNC:
    input: "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb"
    output:
        ECNC_file = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/ECNC/{genome}_RecBlastOutput.bed.rbb.ecnc",
        dropped_lines = "output/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/ECNC/dropped/{genome}_RecBlastOutput.bed.rbb.ecnc.dropped"
    log: "logs/ECNC/{genome}/{query}-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/{genome}_RecBlastOutput.bed.rbb.ecnc.log"
    conda: "../envs/conda_recBlat.yaml"
    script: "../code/ecnc_sm.py"

#
# rule GetGeneSequence:
#     input:
#         genome = "data/genome/{genome}.fa",
#         BED_file = "output/geneBedFiles/{genome}_{gene}.bed"
#     output:
#         "output/geneSequences/{gene}_{genome}_{upstream}kbUp_{downstream}kbDown.fa"
#     shell: "bedtools | getFastaFromBed -name -s -fi <input FASTA>"
