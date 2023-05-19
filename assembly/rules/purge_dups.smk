rule pd_split_assembly:
    input: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa"
    output: "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/split.fa"
    conda: "../envs/HiFiAssembly.yml"
    shell: "split_fa {input} > {output}"


rule pd_merge_CCS:
    input: get_hifiasm_inputs
    output: "output/HiFi-adapterFiltered-merged/{settings}/{species}.ccs.filt.fastq.gz"
    conda: "../envs/HiFiAssembly.yml"
    shell: "cat {input} > {output}"

rule pd_map_CCS_to_assembly:
    input:
        genome = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/split.fa",
        fastq = "output/HiFi-adapterFiltered-merged/{ccs_opts}/{species}.ccs.filt.fastq.gz"
    output: "output/minimap2/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.HiFi.paf"
    conda: "../envs/HiFiAssembly.yml"
    threads: 32
    shell: "minimap2 -t {threads} -I 300G -x map-hifi {input.genome} {input.fastq} > {output}"

rule pd_map_assembly_to_assembly:
    input:
        genome = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/split.fa"
    output: "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.genome.paf"
    conda: "../envs/HiFiAssembly.yml"
    threads: 32
    shell: "minimap2 -t {threads} -x asm5 -I 300G -DP {input.genome} {input.genome} > {output}"


rule pd_haploid_diploid_coverage:
    input: "output/minimap2/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.HiFi.paf"
    output: multiext("output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/PB", ".stat" ,".base.cov")
    conda: "../envs/HiFiAssembly.yml"
    params: 
        outdir = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/"
    shell: "pbcstat -O {params.outdir} {input} "


rule pd_calcuts_assembly:
    input: "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/PB.stat"
    output: "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/cutoffs"
    log: "logs/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/calcuts.log"
    params: 
        outdir = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/"
    conda: "../envs/HiFiAssembly.yml"
    shell: "ROOTDIR=$(pwd -P) ; "
           "cd {params.outdir} ; "
           "calcuts $ROOTDIR/{input} > $ROOTDIR/{output} 2> $ROOTDIR/{log} ; "
           "cd - "

rule purge_dups:
    input: 
        cov = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/PB.base.cov",
        cutoffs = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/cutoffs", 
        ava = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.genome.paf"
    output: 
        dups = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/dups.bed",
    params: 
        outdir = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/"
    conda: "../envs/HiFiAssembly.yml"
    shell: "ROOTDIR=$(pwd -P) ; "
           "cd {params.outdir} ; "
           "purge_dups -2 -c $ROOTDIR/{input.cov} -T $ROOTDIR/{input.cutoffs} $ROOTDIR/{input.ava} > dups.bed ; "
           "cd -"

rule pd_getseq:
    input: 
        dups = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/dups.bed",
        og = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa",
    output: 
        fa = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}-purgedDups.fa",        
        ln = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}-purgedDups.fa",        
    params: 
        outdir = "output/purge_dups/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/"
    conda: "../envs/HiFiAssembly.yml"
    shell: "ROOTDIR=$(pwd -P) ; "
           "cd {params.outdir} ; "
           "get_seqs -e -p {wildcards.species}.{wildcards.genometype}.{wildcards.hic} $ROOTDIR/{input.dups} $ROOTDIR/{input.og} ; "
           "ln -sf $ROOTDIR/{output.fa} $ROOTDIR/{output.ln} ; "
           "cd -"


#Map reads to assembly and assembly to self
#./purge_dups/src/split_fa input.fasta > split.fasta

#./minimap2 -I 200G -t 24 -xasm5 -DP split.fasta split.fasta > split.genome.paf

#./minimap2 -I 200G -x map-pb -t 24 split.fasta reads.ccs.fastq.gz > reads.paf

#Calculate haploid/diploid coverage threshold and remove haplotype duplicates from assembly
#./purge_dups/src/pbcstat -O coverage reads.paf

#./purge_dups/src/calcuts PB.stat > cutoffs

#./purge_dups/src/purge_dups -2 -c PB.base.cov -T cutoffs split.genome.paf > dups.bed

#./purge_dups/src/get_seqs -e -p asm_mTadBra.purged dups.bed input.fasta



