rule remove_problem_genes_updated:
    version: 2
    input:
        fa = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/predict_results/{species}.proteins.fa", 
        problems = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/predict_results/{species}.models-need-fixing.txt"
    output: 
        "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/predict_results/{species}.proteins.noProblemGenes.fa"
    run: 
        from Bio import SeqIO
        with open(input.problems) as errors:
            remove_me = [i.split('\t')[0] for i in errors if not(i.startswith('#'))]
        keep = []
        for record in SeqIO.parse(input.fa, "fasta"):
            gene = record.description.split(" ")[1]
            if gene in remove_me:
                continue
            else:
                keep.append(record)
        with open(output[0], "w") as outhandle:
            SeqIO.write(keep, outhandle, "fasta")


rule split_prot_updated:
    version: 2
    input:
        "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/predict_results/{species}.proteins.noProblemGenes.fa"
    output:
        expand("output/funannotate/TOGA_Prot/{{genome}}_TOGA_Prot/annotate_misc/proteins_split/{{species}}.proteins.part_{part}.fa", part=[str(i) + "fa" for i in [j for j in range(1,21)]])
    params:
        parts = 20,
        dir = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/proteins_split"
    conda: 
        "../envs/liftoff.yaml"
    threads: 32
    shell:
        "seqkit split2 "
        " -p {params.parts} "
        " --out-dir {params.dir} "
        " -t {threads} "
        " {input}; "
        "rename '_00' '_' {output.dir}/*.fa; "
        "rename '_0' '_' {output.dir}/*.fa"
        

#rule interproscan_manual_split:
#    input: 
#        proteins = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/proteins_split/{species}.proteins.part_{part}.fa",
#    output: 
#        xml = temp("output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/{species}.proteins.part_{part}.fa.xml")
#    params:
#        outdir = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/"
#    threads: 32
#    shell: 
#        "module load java && "
#        "module load gcc && "
#        " $(pwd -P)/code/interproscan-5.60-92.0/interproscan.sh "
#        " --cpu {threads} "
#        " -i $(pwd -P)/{input.proteins} "
#        " XML " 
#        " -d $(pwd -P)/{params.outdir} "
#        " -goterms -pa; "
#        " module unload gcc;"
#        " module unload java"


rule interproscan_manual:
    input: 
        proteins = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/predict_results/{species}.proteins.noProblemGenes.fa"
    output: 
        xml = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/{species}.proteins.noProblemGenes.fa.xml"
    params:
        outdir = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/"
    threads: 32
    shell: 
        "module load java && "
        "module load gcc && "
        " $(pwd -P)/code/interproscan-5.60-92.0/interproscan.sh "
        " --cpu {threads} "
        " -i $(pwd -P)/{input.proteins} "
        " XML " 
        " -d $(pwd -P)/{params.outdir} "
        " -goterms -pa; "
        " module unload java; "
        " module unload gcc; "


rule interproscan_split_join:
    input: 
        expand("output/funannotate/TOGA_Prot/{{genome}}_TOGA_Prot/annotate_misc/{{species}}.proteins.part_{part}.fa.xml", part= [str(i) for i in range(1,21)])
    output: 
        xml = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/{species}.proteins.noProblemGenes.fa.xml"
    params:
        outdir = "output/funannotate/TOGA_Prot/{genome}_TOGA_Prot/annotate_misc/"
    threads: 1
    shell: 
        "echo '<protein-matches>' > {output}; "
        "cat {input} | "
        " sed '/protein-matches/d' >> "
        " {output}; "
        "echo '</protein-matches>' >> {output}"
