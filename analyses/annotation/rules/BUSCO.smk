rule BUSCO_prot_mammalia:
    input: "output/funannotate/{genome}/update_results/{species}.proteins.fa"
    output: directory("output/BUSCO/{genome}/{species}_proteome_mammalia/")
    threads: 32
    log:
        stdout = "logs/BUSCO_prot_mammalia/{genome}/{species}_BUSCO.stdout.txt",
        stderr = "logs/BUSCO_prot_mammalia/{genome}/{species}_BUSCO.stderr.txt"
    shell:
        "singularity run busco_v5.4.3_cv1.sif busco -f "
        " --in {input} "
        " --out {output} "
        " -l mammalia_odb10 "
        " --mode proteins "
        " -c {threads} "
        " 1> {log.stdout} "
        " 2> {log.stderr}"


rule BUSCO_gff_evidence_eutheria:
    input: "data/GFF_evidences_protein/{type}/{file}.fa"
    output: directory("output/BUSCO/GFF_evidence/{type}/{file}_eutheria/")
    threads: 32
    log:
        stdout = "logs/BUSCO/GFF_evidence/{type}/{file}_eutheria.stdout.txt",
        stderr = "logs/BUSCO/GFF_evidence/{type}/{file}_eutheria.stderr.txt"
    shell:
        "singularity run busco_v5.4.3_cv1.sif busco -f "
        " --in {input} "
        " --out {output} "
        " -l mammalia_odb10 "
        " --mode proteins "
        " -c {threads} "
        " 1> {log.stdout} "
        " 2> {log.stderr}"
	
	
