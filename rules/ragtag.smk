rule ragtag_haps:
    version: 0.1
    input: 
        ref = "data/assemblies/{ccs_settings}/{hifi_opts}/{species}.p_ctg.hic_mapq0_noContigEC_minQ10_scaffolds_final_manualCuration.chr_M.fa",
        hap = "data/assemblies/{ccs_settings}/{hifi_opts}/{species}.{hap}.p_ctg.hic.chr_M.fa"
    threads: 20
    output: directory("output/ragtag_hap/{ccs_settings}/{hifi_opts}/{species}_{hap}")
    shell: "ragtag_scaffold.py "
           "{input.ref} "
           "{input.hap} "
           "--aligner nucmer "
           "-u "
           "-o {output}"

rule link_ragtag:
    input: "output/ragtag_hap/{ccs_settings}/{hifi_opts}/{species}_{hap}/ragtag.scaffold.fasta"
    output: "data/assemblies/{ccs_settings}/{hifi_opts}/{species}.{hap}.p_ctg.hic_mapq0_noContigEC_minQ10_scaffolds_final_manualCuration.chr_M.ragtag_scaffolded.fa"
    shell: "PROJDIR=$(pwd -P); "
           "ln -s $PROJDIR/{input} $PROJDIR/{output}"

