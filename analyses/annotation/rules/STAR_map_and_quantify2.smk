# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
mapping to a STAR index

Expects a json config file with the following structure, 

{
  "input_samples": {
    "mouse-a-brain_brain": [
      "SRR594393_1.fastq.gz", 
      "SRR594393_2.fastq.gz"
    ], 
    "mouse-c-testes_testes": [
      "SRR594418_1.fastq.gz", 
      "SRR594418_2.fastq.gz"
    ]
  }, 
  "map_opts": "", 
  "species": "mus_musculus"
}

"""

__author__ = "Peter Sudmant"
__license__ = "MIT"

#import pysam 
import os
import tempfile
import pandas as pd
import pdb

configfile: "config_trimmed.json"

def get_rsem_inputs_w_genenames(wildcards):
    inputs = []
    for sample in [k for k in config['input_samples'].keys()]: 
        inputs.append("output/rsem/{sample}.genes.w_names.results".format(sample=sample))
    return inputs        

def get_rsem_inputs(wildcards):
    inputs = []
    for sample in [k for k in config['input_samples'].keys()]: 
        inputs.append("output/rsem/{sample}.isoforms.results".format(sample=sample))
    return inputs        

def get_map_inputs(wildcards):
    inputs = []
    for sample in [k for k in config['input_samples'].keys()]: 
        inputs.append("output/mapping/{sample}/Aligned.sortedByCoord.out.bam.bai".format(sample=sample))
    return inputs        

def get_inputs(wildcards):
    inputs = []
    inputs += get_map_inputs(wildcards)
    #inputs += get_rsem_inputs(wildcards)
    if config['rsem_skip'] != 1: 
        inputs += ["output/rsem_summary/genes.results"]
        inputs += ["output/rsem_summary/isoforms.results"]
    return inputs        

rule STAR_all:
    input:
        get_inputs
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem-per-cpu 1000 "
                                     "--mem 4000 "
                                     "--time 0-0:05:00 "
                                     "-J rule_all "
                                     "-o logs/rule_all_%j.logs "
                                     "-p defq "

rule RSEM_TMM:
    output:
        temp("output/rsem_summary/tmp.genes.tmm.results"),
        "output/rsem_summary/genes.tmm.results"
    input:
        "output/rsem_summary/genes.tpm.results"
    run:
        shell("Rscript /home/psudmant/code/seqlib/seqlib/expression/TMM_normalize.R --fn_input {input[0]} --fn_output {output[0]}")
        t = pd.read_csv(output[0], header=0, sep="\t")
        cols = list(t.columns)
        cols[0] = "gene_name"
        cols = [c[0] == "X" and c[1:] or c for c in cols]
        t.columns = cols
        t.to_csv(output[1], sep="\t", index=False)


def get_gene_id_table(fn_gtf):
    gtf_cols=["contig",
              "source", 
              "feature", 
              "start", 
              "end", 
              "score", 
              "strand", 
              "frame", 
              "attribute"]
    t = pd.read_csv(fn_gtf, sep="\t",names = gtf_cols)
    t.attribute = t.attribute.str.replace('"','')
    IDs_to_names = pd.DataFrame(t.attribute.str.split(";*\s*").tolist())
    row_0 = IDs_to_names.iloc[0].values.tolist()
    #pdb.set_trace()
    idx_tid, idx_gid, idx_gname = [row_0.index(id)+1 for id in ["transcript_id","gene_id","gene_name"]]
    IDs_to_names = IDs_to_names[[idx_tid, idx_gid, idx_gname]]
    IDs_to_names.columns = ['transcript_id', 'gene_id', 'gene_name']
    IDs_to_names = IDs_to_names[['gene_id', 'gene_name']].drop_duplicates()
    
    return IDs_to_names

rule RSEM_summary:
    output:
        "output/rsem_summary/{type}.results"
    input:
        lambda wildcards: ["output/rsem/{sample}.{type}.results".format(sample=sample,type=wildcards.type) for sample in list(config['input_samples'].keys())]
    threads:
        1
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 6000 "
                                     "--time 0-8:00:00 "
                                     "-J make_RSEM_summary "
                                     "-o logs/make_RSEM_summary_%j.logs "
                                     "-p defq "
    run:
        gene_ids = get_gene_id_table(config['fn_gtf'])
        
        tables = []
        for fn in input:
            splt = ".{type}.results".format(type=wildcards.type)
            sample = fn.split("/")[-1].split(splt)[0]
            t = pd.read_csv(fn, header=0, sep="\t")
            t['sample'] = sample
            tables.append(t)
        t = pd.concat(tables)
        t.gene_id = t.gene_id.astype("str") #if all #s it casts this as a double or something
        t = pd.merge(t, gene_ids, left_on='gene_id', right_on='gene_id')
        
        t.to_csv(output[0], index=False, sep="\t") 
    
rule RSEM:
    output:
        "output/rsem/{sample}.isoforms.results",
        "output/rsem/{sample}.genes.results"
    input:
        "output/mapping/{sample}/Aligned.toTranscriptome.out.bam"
    threads: 26
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 100000 "
                                     "--time 0-8:00:00 "
                                     "-J RSEM_{sample} "
                                     "-o logs/RSEM_{sample}_%j.logs "
                                     "-p defq ".format(sample=wildcards.sample)
    run:
        tmp_dir = "/dev/shm/RSEM_{sample}_tmp/".format(sample=wildcards.sample)

        pe = "--paired-end "
        if config['type'] == "se":
            pe = " "
        #print(config['reference_by_species'])
        shell("rsem-calculate-expression  "
               "-p {{threads}} "
               "--no-bam-output "
               "--temporary-folder {tmp_dir} "
               "{rsem_opts} "
               "--bam {pe} {fn_bam} "
               "{rsem_path}/{annot}/{species}/{ref}/{ref} "
               "output/rsem/{sample} ".format(#RSEM_threads=24,
                                       rsem_path = config['rsem_path'],
                                       annot = config['annotation'],
                                       species = config['species'],
                                       ref = config['reference_by_species'][config['species']],
                                       rsem_opts = config['rsem_opts'],
                                       pe=pe,
                                       tmp_dir = tmp_dir,
                                       sample=wildcards.sample,
                                       fn_bam = input[0]))

rule index_bams:
    input:
        "output/mapping/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "output/mapping/{sample}/Aligned.sortedByCoord.out.bam.bai"
    threads: 13
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 30000 "
                                     "--time 0-6:00:00 "
                                     "-J index_{sample} "
                                     "-o logs/index_{sample}_%j.logs "
                                     "-p defq "
                                     "".format(sample=wildcards.sample)
    shell:
        "samtools index -@ {threads} {input[0]}"

def get_map_fastq_inputs(wildcards):
    inputs = []
    
    for fastq in config['input_samples'][wildcards.sample]:
        inputs.append("{sample_path}/"
                       "{sample_fastq}"
                       "".format(sample_path=config['input_sample_path'],
                                 sample=wildcards.sample,
                                 sample_fastq=fastq))
    #print(inputs)
    return inputs

rule map:
    input:
        get_map_fastq_inputs,
        "{index_path}/SA".format(index_path=config['STAR_index'])
    output:
        "output/mapping/{sample}/Aligned.sortedByCoord.out.bam",
        "output/mapping/{sample}/Aligned.toTranscriptome.out.bam"
    shadow: "shallow"
    threads: 10
    run:
        n_files = len(input)-1
        all_files = sorted(input[:n_files])

        if config['type'] == "se":
            read1 = ",".join(all_files)
            read2= ""
        else:
            assert n_files % 2 == 0
            read1s, read2s = [], []
            
            for i in range(int(n_files/2)):
                read1s.append(all_files[i*2])
                read2s.append(all_files[i*2+1])
            read1 = ",".join(read1s)     
            read2 = ",".join(read2s)     
        # REALLY important that you set 'dir=' here to be your scratch folder        
        with tempfile.TemporaryDirectory(prefix="STAR_%s_"%(wildcards.sample), dir="/global/scratch2/mvazquez/") as tmp_dir:
            shell("STAR runThreadN {STAR_threads} "
                    "--genomeDir {index_path} "
                    "--readFilesIn {read1} {read2} "
                    "--readFilesCommand zcat "
                    "{map_opts} "
                    "--outTmpDir {tmp_dir}/tmp "
                    "--outFileNamePrefix output/mapping/{sample}/ "
                    "--outSAMtype BAM SortedByCoordinate  "
                    "--quantMode GeneCounts TranscriptomeSAM "
                    "--limitBAMsortRAM 60000000000 "
                    "--sjdbGTFfile {fn_gtf} "
                    " ".format(index_path=config['STAR_index'],
                                          STAR_threads='{threads}',
                                          sample=wildcards.sample,
                                          fn_gtf=config['fn_gtf'],
                                          tmp_dir = tmp_dir,
                                          map_opts = config['map_opts'],
                                          read1=read1,
                                          read2=read2))
        """cleanup the tempdir"""
        #tmp_dir.cleanup()

rule clean:
    shell:
        "rm -rf output/mapping; rm logs/*; rm -rf output/rsem"

