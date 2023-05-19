## Note: this fails when an assembly - always a haplotig - lacks a mitogenome. In these cases, I manually create chr_M.fa by copying it from the consensus chr_M and do the last steps. 

rule one_mitochondria: 
    input: 
        fa = "data/assemblies/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}.fa",
        fai = "data/assemblies/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}.fa.fai",
        mitofa = "output/mitoHiFi/fromContig/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/final_mitogenome.fasta",
        mitocontigs = "output/mitoHiFi/fromContig/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/contigs_stats.tsv"
    output:
        fa = "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.chr_M.fa",
        nonMitofa = "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.noMito.fa",
        mitofa = "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.mitocontigs.fa",
        mitocontigs = "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/mitocontigs.txt",
        nonMitocontigs = "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/non-mitocontigs.txt",
        chr_M = "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/chr_M.fa",
    params:
        mitocontigs = "output/mitoHiFi/fromContig/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/potential_con",
        outdir = "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/"
    shell: "mkdir -p {params.outdir} ; "
           "cat {input.mitocontigs} | cut -f1 | grep -v 'contig_id' | grep -v 'final_mitogenome' > {output.mitocontigs} || echo 'No other mitocontigs, final only!'; "
           "awk '{{print $1}}' {input.fai} | grep -v -f {output.mitocontigs} > {output.nonMitocontigs} || echo 'failed line 3'; "
           "samtools faidx -o {output.nonMitofa} {input.fa} -r {output.nonMitocontigs} || echo 'failed line 4'; "
           "cat {input.mitofa} | sed 's/^>.*$/>chr_M/' > {output.chr_M} || echo 'failed line 5'; "
           "cat {output.nonMitofa} {output.chr_M} > {output.fa} || echo 'failed line 6' ; "
           "samtools faidx -o {output.mitofa} {input.fa} -r {output.mitocontigs} || echo 'no other mitocontigs, final only!'"


rule link_oneMito_fa:
    input: "output/genomes_curated/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.chr_M.fa"
    output: "data/assemblies/{ccs_opts}/{hifi_opts}/{species}.{genometype}.{hic}.chr_M.fa"
    shell: "ROOTDIR=$(pwd -P); "
           "ln -sf $ROOTDIR/{input} $ROOTDIR/{output}"


rule wfmash_align_genomes_to_ref:
    input: 
        genome1 = "data/assemblies/{ccs_opts}/{hifi_opts}/{genome1}.fa",
        genome2 = "data/assemblies/{ccs_opts}/{hifi_opts}/{genome2}.fa"
    output: "output/wfmash/aligned_ref/{ccs_opts}/{hifi_opts}/t_{genome1}-q_{genome2}.approx.fa"
    params:
        kmer=21
    threads: 40
    shell: "wfmash {genome1} {genome2} -t{threads} -k{params.kmer} -m > {output}"
