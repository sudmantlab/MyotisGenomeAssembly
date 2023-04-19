rule liftoff:
    input: 
        refGenome="data/genomes/{refGenome}.fa",
        gff="data/GFF/{refGenome}_{target}.gff",
        queryGenome="data/genomes/{genome}.fa",
    output: 
        placed="output/liftoff/{refGenome}-{genome}-{target}.gff",
        unplaced="output/liftoff/{refGenome}-{genome}-{target}.unplaced",
        tmp=temp(directory("output/liftoff/tmp_{refGenome}_{genome}_{target}"))
    threads: 40
    conda: "../envs/liftoff.yaml"
    shell: "liftoff -g {input.gff} -u {output.unplaced} -p 40 -exclude_partial -cds -polish -dir {output.tmp} -o {output.placed} {input.queryGenome} {input.refGenome}"	

rule liftoff_getfasta:
    input: 
        fi = "data/genomes/{genome}.fa",
        bed = "output/liftoff/{refGenome}-{genome}-{target}.gff"
    output: "output/liftoff_fasta/{genome}-{refGenome}-{target}.fa"
    conda: "../envs/liftoff.yaml"
    shell: "bedtools getfasta "
           "-fi {input.fi} "
           "-name+ -split "
           "-bed {input.bed} "
           "-fo {output}"
