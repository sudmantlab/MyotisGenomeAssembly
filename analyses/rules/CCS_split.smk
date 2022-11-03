import os


rule ccs_consensus_noopts_chunk:
    input: "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{species}/defaults-minPasses3_minRQ0.99/{pacbio1}/{pacbio2}/{id}.ccs.{chunk}_40.bam",
        report="output/HiFi-CCS/{species}/defaults-minPasses3_minRQ0.99/{pacbio1}/{pacbio2}/{id}.ccs_report.{chunk}_40.txt"
    log: "logs/HiFi-CCS/{species}/defaults-minPasses3_minRQ0.99/{pacbio1}/{pacbio2}/{id}.{chunk}.logs"
    params:
        loglevel = "DEBUG"
    wildcard_constraints:
        chunk="[0-9]?[0-9]",
        species="[A-Z]_[a-z]+"
    threads: 4
    conda: "../envs/HiFiAssembly.yml"
    shell: "ccs -j {threads} --min-passes {wildcards.minPasses} --min-rq {wildcards.minRQ} --report-file {output.report} --log-level {params.loglevel} --log-file {log} --chunk {wildcards.chunk}/40 {input} {output.ccs}"


rule ccs_consensus_minPasses_minRQ_chunk:
    input: "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{species}/minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs.{chunk}_40.bam",
        report="output/HiFi-CCS/{species}/minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs_report.{chunk}_40.txt"
    log: "logs/HiFi-CCS/{species}/minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.{chunk}.logs"
    params:
        loglevel = "DEBUG"
    wildcard_constraints:
        chunk="[0-9]?[0-9]",
        species="[A-Z]_[a-z]+"
    threads: 4
    conda: "../envs/HiFiAssembly.yml"
    shell: "ccs -j {threads} --min-passes {wildcards.minPasses} --min-rq {wildcards.minRQ} --report-file {output.report} --log-level {params.loglevel} --log-file {log} --chunk {wildcards.chunk}/40 {input} {output.ccs}"

rule ccs_chunk_merge:
    input: expand("output/HiFi-CCS/{species}/{isDefault}minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs.{chunk}_40.bam", chunk=range(1,41), allow_missing=True)
    output: "output/HiFi-CCS/{species}/{isDefault}minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs.bam"
    wildcard_constraints:
        isDefault="(defaults-)?"
        chunk="\d{1,2}",
        species="[A-Z]_[a-z]+"
    threads: 32
    shell: "samtools merge -@{threads} {output} {input}"
    

rule HiFiAdapterFilt:
    input: "output/HiFi-CCS/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.bam"
    output: 
        blocklist = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.blocklist", 
        contaminant = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.contaminant.blastout", 
        fasta = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.fasta", 
        fastq = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.fastq",
        filtered = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq"
    log: "logs/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.log"
    params:
        outDir = lambda wildcards, output: os.path.dirname(output[0]),
        inPrefix = lambda wildcards, input: input[0].rstrip(".bam")
    conda: "../envs/HiFiAssembly.yml"
    threads: 32
    shell: "pbadapterfilt.sh -b {params.inPrefix} -t {threads} -o {params.outDir} &> {log}"


rule hifiasm_noopts:
    input: "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq"
    output: 
        r_utg = "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm.r_utg.gfa",
        p_utg = "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm.p_utg.gfa",
        p_ctg = "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm.p_ctg.gfa",
        a_ctg = "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm.a_ctg.gfa"
    params:
        prefix = "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm"
    threads: 32
    log: "output/hifiasm/{species}/{pacbio1}/{pacbio2}/{id}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "hifiasm -o {params.prefix} -t {threads} {input} > {log} 2>&1"
