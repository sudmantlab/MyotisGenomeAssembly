import os
import tempfile
import pandas as pd
import pdb

#configfile: "config.yaml"


#def get_STAR_all(wildcards):
#    inputs = []
#    out_path= "output/mapping/{species}/{genome}-GFF/{sample}/Aligned.sortedByCoord.out.bam.bai"
#    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
#    for species in config['species']:
#        #print(species)
#        for ref in config['reference_by_species'][species]:
#            #print(ref)
#            inputs.extend([out_path.format(species = species, genome = ref, sample = r[3]) for r in samples[samples["Species"] == species].itertuples()])
#            #print(inputs)
#    return inputs


#rule STAR_all:
#    input:
#        get_STAR_all


#def get_rsem_inputs_w_genenames(wildcards):
#    '''Cavet emptor, this is Peter's original code.'''    
#    inputs = []
#    for sample in [k for k in config['input_samples'].keys()]:
#        inputs.append("output/rsem/{sample}.genes.w_names.results".format(sample=sample))
#    return inputs

#def get_rsem_inputs(wildcards):
#    '''Cavet emptor, this is Peter's original code.'''    
#    inputs = []
#    for sample in [k for k in config['input_samples'].keys()]:
#        inputs.append("output/rsem/{sample}.isoforms.results".format(sample=sample))
#    return inputs


#rule RSEM_TMM:
#    '''Cavet emptor, this is Peter's original code.'''    
#    output:
#        temp("output/rsem_summary/tmp.genes.tmm.results"),
#        "output/rsem_summary/genes.tmm.results"
#    input:
#        "output/rsem_summary/genes.tpm.results"
#    run:
#        shell("Rscript /home/psudmant/code/seqlib/seqlib/expression/TMM_normalize.R --fn_input {input[0]} --fn_output {output[0]}")
#        t = pd.read_csv(output[0], header=0, sep="\t")
#        cols = list(t.columns)
#        cols[0] = "gene_name"
#        cols = [c[0] == "X" and c[1:] or c for c in cols]
#        t.columns = cols
#        t.to_csv(output[1], sep="\t", index=False)


def get_gene_id_table(fn_gtf):
    '''Cavet emptor, this is Peter's original code.'''
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
    idx_tid, idx_gid, idx_gname = [row_0.index(id)+1 for id in 
                                   ["transcript_id","gene_id","gene_name"]]
    IDs_to_names = IDs_to_names[[idx_tid, idx_gid, idx_gname]]
    IDs_to_names.columns = ['transcript_id', 'gene_id', 'gene_name']
    IDs_to_names = IDs_to_names[['gene_id', 'gene_name']].drop_duplicates()

    return IDs_to_names

#rule RSEM_summary:
#    '''Caveat emptor, this is Peter's original code.'''
#    output:
#        "output/rsem_summary/{type}.results"
#    input:
#        lambda wildcards: 
#            ["output/rsem/{species}/{genome}-GFF/{sample}.{type}.results".format(sample=sample,type=wildcards.type) 
#             for sample in list(config['input_samples'].keys())]
#    threads: 1
#    run:
#        gene_ids = get_gene_id_table(config['fn_gtf'])
#
#        tables = []
#        for fn in input:
#            splt = ".{type}.results".format(type=wildcards.type)
#            sample = fn.split("/")[-1].split(splt)[0]
#            t = pd.read_csv(fn, header=0, sep="\t")
#            t['sample'] = sample
#            tables.append(t)
#        t = pd.concat(tables)
#        t.gene_id = t.gene_id.astype("str") #if all #s it casts this as a double or something
#        t = pd.merge(t, gene_ids, left_on='gene_id', right_on='gene_id')
#
#        t.to_csv(output[0], index=False, sep="\t")


rule RSEM:
    output:
        "output/rsem/{species}/{ref}-{type}-{id}/{sample}.isoforms.results",
        "output/rsem/{species}/{ref}-{type}-{id}/{sample}.genes.results",
        directory("output/rsem/{species}/{ref}-{type}-{id}/{sample}.stat")
    input:
        bam = "output/mapping/{species}/{ref}-{type}-{id}/{sample}/Aligned.toTranscriptome.out.bam",
        gff = "data/GFF_evidences/{type}/{ref}_{id}.gff3",
        rsem = "output/RSEM_indexes/{species}/{ref}-{type}-{id}/{ref}.ngvec"
    threads: 56
    params:
        rsem_opts = '' #config['rsem_opts']
    shadow: 'full'
    conda: '../envs/STAR-RSEM-EBSeq.yaml'
    shell: "rsem-calculate-expression  "
           "-p {threads} "
           "--no-bam-output "
           "{params.rsem_opts} "
           "--bam "
           "--paired-end "
           "{input.bam} "
           "output/RSEM_indexes/{wildcards.species}/{wildcards.ref}-{wildcards.type}-{wildcards.id}/{wildcards.ref} "
           "output/rsem/{wildcards.species}/{wildcards.ref}-{wildcards.type}-{wildcards.id}/{wildcards.sample} "

rule index_bams_withGFF:
    input:
        "output/mapping/{species}/{ref}-{type}-{id}/{sample}/Aligned.{toType}.out.bam"
    output:
        "output/mapping/{species}/{ref}-{type}-{id}/{sample}/Aligned.{toType}.out.bam.bai"
    threads: 56
    conda: '../envs/STAR-RSEM-EBSeq.yaml'
    shell:
        "samtools index -@ {threads} {input[0]}"


def get_map_fastq_inputs_withGFF(wildcards):
    path_trimmed = "data/RNA-seq_trimmed/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    import pandas as pd
    from pathlib import Path
    samples = pd.read_table("rna_pepsamples.tsv", index_col= False)
    samples = samples[samples.Sample == wildcards.sample]
    input_dict = {"left": [], 
                  "right": [], 
                  "index": "output/STAR_indexes/{species}/{ref}-{type}-{id}/SA",
                  "gff": "data/GFF_evidences/{type}/{ref}_{id}.gff3"}
    samples_grouped = samples.groupby(samples.Sample)
    for sample in set(samples["Sample"].tolist()):
        sample_subset = samples_grouped.get_group(sample)
        #print(sample_subset)
        if len(sample_subset) == 0:
            raise Exception("No files available for sample {}".format(sample))
        #print(r[2])
        input_dict["left"].extend([path_trimmed.format(species = r[2], 
                                                       sample_name = r[3], 
                                                       read = r[4]) 
                                   for r in 
                                   sample_subset[sample_subset["Read"] == "R1"].itertuples()
                                   ])
        input_dict["right"].extend([path_trimmed.format(species = r[2], 
                                                        sample_name = r[3], 
                                                        read = r[4]) 
                                    for r in 
                                    sample_subset[sample_subset["Read"] == "R2"].itertuples()
                                    ])
    print(input_dict)
    return input_dict


rule map_pe_withGFF:
    input: unpack(get_map_fastq_inputs_withGFF)
    output:
        "output/mapping/{species}/{ref}-{type}-{id}/{sample}/Aligned.sortedByCoord.out.bam",
        "output/mapping/{species}/{ref}-{type}-{id}/{sample}/Aligned.toTranscriptome.out.bam",
        "output/mapping/{species}/{ref}-{type}-{id}/{sample}/Log.final.out"
    shadow: "shallow"
    log: "logs/mapping/{species}/{ref}-{type}-{id}/{sample}/star.log"
    params:
        map_opts = '' #config['map_opts']
    conda: "../envs/STAR-RSEM-EBSeq.yaml"
    threads: 56
    shell:
           "STAR runThreadN {threads} "
           "--genomeDir output/STAR_indexes/{wildcards.species}/{wildcards.ref}-{wildcards.type}-{wildcards.id} "
           "--readFilesIn {input.left} {input.right} "
           "--readFilesCommand zcat "
           "--sjdbGTFfile {input.gff} "
           "--sjdbGTFtagExonParentTranscript Parent "
           "{params.map_opts} "
           "--outFileNamePrefix output/mapping/{wildcards.species}/{wildcards.ref}-{wildcards.type}-{wildcards.id}/{wildcards.sample}/ "
           "--outSAMtype BAM SortedByCoordinate  "
           "--twopassMode Basic "
           "--quantMode GeneCounts TranscriptomeSAM "
           "--limitBAMsortRAM 60000000000 "
           "> {log} 2>&1"


