rule pretext_map:
    input: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.sorted.bam"
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pretext"
    conda: "../envs/RapidCuration.yaml"
    threads: 2
    shell: "samtools view -h {input} | "
           " awk '{{print gensub(/scaff_([0123456789]+)_contig_tig([0123456789]+)[^\t]+(\s)/,\"\\1_\\2\\3\",\"g\")}}' | "
           " PretextMap -o {output} --sortby 'length' --sortorder 'descend' --mapq 10 --highRes"
