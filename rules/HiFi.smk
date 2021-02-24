rule ccs_consensus_noopts:
    input: "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{species}/defaults/{pacbio1}/{pacbio2}/{id}.ccs.bam",
        report="output/HiFi-CCS/{species}/defaults/{pacbio1}/{pacbio2}/{id}.ccs_report.txt"
    log: "logs/HiFi-CCS/{species}/defaults/{pacbio1}/{pacbio2}/{id}.logs"
    threads: 32
    conda: "../envs/HiFiAssembly.yml"
    shell: "ccs -j {threads} --report-file {output.report} --log-file {log} {input} {output.ccs}"


rule ccs_consensus_minPasses_minRQ:
    input: "data/PacBio-HiFi/{species}/data2/pb/{pacbio1}/{pacbio2}/{id}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{species}/minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs.bam",
        report="output/HiFi-CCS/{species}/minPasses{minPasses}_minRQ{minRQ}/{pacbio1}/{pacbio2}/{id}.ccs_report.txt"
    log: "logs/HiFi-CCS/minPasses{minPasses}_minRQ{minRQ}/{species}/{pacbio1}/{pacbio2}/{id}.logs"
    threads: 32
    conda: "../envs/HiFiAssembly.yml"
    shell: "ccs -j {threads} --min-passes {wildcards.minPasses} --min-rq {wildcards.minRQ} --report-file {output.report} --log-file {log} {input} {output.ccs}"


rule bam2fastq:
    input: "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs.bam"
    output: "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs.fastq.gz"
    params:
        prefix = "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs"
    conda: "../envs/HiFiAssembly.yml"
    shell: "bam2fastq -o {params.prefix} {input}"


rule hifiasm_noopts:
    input: "output/HiFi-CCS/{species}/{pacbio1}/{pacbio2}/{id}.ccs.fastq.gz"
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
