rule longQC_bam:
    input: "output/HiFi-CCS/{species}/{settings}/{pb1}/{pb2}/{id}.ccs.fastq"
    output: "output/longQC/{species}/{settings}/{pb1}/{pb2}/{id}/{id}_web_summary.html"
    log: "logs/longQC/{species}/{settings}/{pb1}/{pb2}/{id}.log"
    params:
        in_file = lambda wildcards: "/input/{species}/{settings}/{pb1}/{pb2}/{id}.ccs.fastq".format(species=wildcards.species,
                                                                                                  settings=wildcards.settings, 
                                                                                                  pb1=wildcards.pb1,
                                                                                                  pb2=wildcards.pb2,
                                                                                                  id=wildcards.id),
        out_dir = lambda wildcards: "/output/{species}/{settings}/{pb1}/{pb2}/{id}/".format(species=wildcards.species,
                                                                                           settings=wildcards.settings, 
                                                                                           pb1=wildcards.pb1,
                                                                                           pb2=wildcards.pb2,
                                                                                           id=wildcards.id),
        out_dir2 = lambda wildcards: "output/longQC/{species}/{settings}/{pb1}/{pb2}/{id}/".format(species=wildcards.species,
                                                                                                   settings=wildcards.settings, 
                                                                                                   pb1=wildcards.pb1,
                                                                                                   pb2=wildcards.pb2,
                                                                                                   id=wildcards.id)

    threads: 32
    #singularity: "docker://docmanny/longqc:113f60c_fixed"
    # shell: "sampleqc -s {wildcards.id}_ -x pb-hifi -o /output/{params.out_dir} -p {threads} /input/{params.in_file}"
    shell: """
        mkdir -p {params.out_dir2}
        rmdir {params.out_dir2}
        singularity run -B /global/scratch2 -B $(pwd -P) -B output/HiFi-CCS/:/input -B output/longQC:/output docker://docmanny/longqc:113f60c_fixed sampleqc -x pb-hifi -o {params.out_dir} -p {threads} -s {wildcards.id}_ {params.in_file}
    """
