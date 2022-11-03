def get_qc_files_before_after_trim(wildcards):
    import pandas as pd
    path_trimmed = "output/qc/fastqc-trimmed/{species}/{sample_name}-trimmed_{read}_fastqc.html"
    path_normal = "output/qc/fastqc/{species}/{sample_name}_{read}_001_fastqc.html"
    file_list = []
    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
    file_list = [path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in samples.itertuples()]
    file_list += [path_normal.format(species = r[2], sample_name = r[3], read = r[4]) for r in samples.itertuples()]
    return file_list

rule multiqc:
    input: get_qc_files_before_after_trim
    output: "output/qc/multiqc/multiqc.html"
    #log: "logs/multiqc/multiqc.html"
    params:
        input_dir = ["output/qc/fastqc", "output/qc/fastqc-trimmed"],
        output_dir = "output/qc/multiqc"
    conda: "../envs/QC.yaml"
    shell: "multiqc --force -o {params.output_dir} -n multiqc.html --profile-runtime  {params.input_dir}"


rule fastqc:
    input:
        "data/RNA-seq/{species}/{sample}.fastq.gz"
    output:
        html="output/qc/fastqc/{species}/{sample}_fastqc.html",
        zip="output/qc/fastqc/{species}/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log: "logs/fastqc/{species}/{sample}.log"
    params:
        outdir = "output/qc/fastqc"
    threads: 4
    conda: "../envs/QC.yaml"
#    shadow: "shallow"
    shell: 
        "fastqc --quiet -t {threads} --outdir {params.outdir}/{wildcards.species} {input}  > {log} 2>&1"

rule fastqc_trimmomatic:
    input:
        "output/trimmed/{species}/{sample}.fastq.gz"
    output:
        html="output/qc/fastqc-trimmed/{species}/{sample}_fastqc.html",
        zip="output/qc/fastqc-trimmed/{species}/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log: "logs/fastqc/{species}/{sample}.log"
    params:
        outdir = "output/qc/fastqc-trimmed"
    threads: 4
    conda: "../envs/QC.yaml"
#    shadow: "shallow"
    shell: 
        "fastqc --quiet -t {threads} --outdir {params.outdir}/{wildcards.species} {input}  > {log} 2>&1"

