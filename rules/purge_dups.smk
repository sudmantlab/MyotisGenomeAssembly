rule split_assembly:
    input: lambda wildcards: "output/hifiasm-fasta/{species}/{settings}/{opt}/{species}.p_ctg.gfa" if wildcards.hic == "bp" else "output/hifiasm-fasta-HiC/{species}/{settings}/{opt}/{species}.asm.hic.p_ctg.gfa"
    output: "output/hifiasm-splitFasta/{species}/"


rule minimap2_CCS_assembly:
    input:
        genome = lambda wildcards: "output/hifiasm/{species}/{settings}/{opt}/{species}.asm.bp.p_ctg.gfa" if wildcards.hic == "bp" else "output/hifiasm-HiC/{species}/{settings}/{opt}/{species}.asm.hic.p_ctg.gfa"
        fastq = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz"
    output: "output/minimap2/{species}/{settings}/{opt}/{pacbio1}/{pacbio2}/{species}.asm.{hic}.p_ctg-{id}.paf"
    conda: "../envs/minimap2.yaml"
    threads: 32
    shell: "minimap2 -t {threads} -xasm20 {input.genome} {input.fastq} > {output}"

def get_minimap2_CCS_assembly:
    file_path = "output/minimap2/{species}/{settings}/{opt}/{pacbio1}/{pacbio2}/{species}.asm.{hic}.p_ctg-{id}.paf"
    species = wildcard.species
    settings = wildcards.settings
    options = wildcards.opt
    hic = wildcards.hic
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == species]
    samples = samples.to_records(index=False)
    input_samples = []
    for sample in samples:
        in_file = file_path.format(species=s[0], settings = settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3], hic = hic, opt=options)
        input_samples.append(in_file)
    return input_samples

rule pbcstat_CCS_assembly:
    input: get_minimap2_CCS_assembly
    output: multiext("output/pbcstat/CCS-assembly/{species}/{settings}/{opt}/{hic}/PB", ".stat" ,".base.cov")
    params: 
        outdir = "output/pbcstat/CCS-assembly/{species}/{settings}/{opt}/{hic}/".format(**wildcards)
    shell: "pbcstat -O {params.outdir} {input} "


rule calcuts_CCS_assembly:
    input: "output/pbcstat/CCS-assembly/{species}/{settings}/{opt}/{hic}/PB.stat"
    output: "output/pbcstat/CCS-assembly/{species}/{settings}/{opt}/{hic}/cutoffs"
    log: "output/pbcstat/CCS-assembly/{species}/{settings}/{opt}/{hic}/calcults.log"
    params: 
        outdir = "output/pbcstat/CCS-assembly/{species}/{settings}/{opt}/{hic}/".format(**wildcards)
    shell: """
           cd {params.outdir}
           calcults PB.stat > cutoffs 2> {log}
           cd -
           """


#for i in $pb_list
#do
#	minimap2 -xasm20 $pri_asm $i | gzip -c - > $i.paf.gz
#	done
#	bin/pbcstat *.paf.gz (produces PB.base.cov and PB.stat files)
#	bin/calcuts PB.stat > cutoffs 2>calcults.log




#Map reads to assembly and assembly to self
./purge_dups/src/split_fa input.fasta > split.fasta

./minimap2 -I 200G -t 24 -xasm5 -DP split.fasta split.fasta > split.genome.paf

./minimap2 -I 200G -x map-pb -t 24 split.fasta reads.ccs.fastq.gz > reads.paf

#Calculate haploid/diploid coverage threshold and remove haplotype duplicates from assembly
./purge_dups/src/pbcstat -O coverage reads.paf

./purge_dups/src/calcuts PB.stat > cutoffs

./purge_dups/src/purge_dups -2 -c PB.base.cov -T cutoffs split.genome.paf > dups.bed

./purge_dups/src/get_seqs -e -p asm_mTadBra.purged dups.bed input.fasta



