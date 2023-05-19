def get_omnic_reads_allSpecies(wildcards):
    import pandas as pd
    path_trimmed = "output/trimmed-hic/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    #input_list = []
    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
    samples = samples[samples.type == "Hi-C"]
    return samples.apply(lambda row: path_trimmed.format(**row), axis=1).tolist()

rule get_readgroupinfo:
    input: get_omnic_reads_allSpecies
    output: 'data/Hi-C/seqinfo.tsv'
    shell: """
        echo -e "fastq\tlibrary\tinstrument_name\trun_id\tflowcell_id\tflowcell_lane" > {output} ; 
        ls {input} | 
          parallel 'echo -en "{{/.}}\t" | 
                      sed "s/.fastq//"; 
                    echo -en "{{/.}}\t" | 
                      sed "s/-trimmed//;s/.fastq//" ; 
                    zcat {{}} | 
                      head -n21 | 
                      grep @ | 
                      cut -d":" -f1,2,3,4 | 
                      sort -u | 
                      sed "s/@//; s/:/\t/g"' | 
          grep _R1 | 
          sed "s/_R1//g" >> {output}
    """


def get_seqinfo(wildcards):
    import pandas as pd
    samples = pd.read_table(config["input_seqinfo"], index_col="fastq")
    #print(samples)
    #print(wildcards.fastq)
    samples = samples.filter(items=[wildcards.fastq], axis=0)
    samples = samples.to_dict('records')[0]
    #print(samples)
    samples["sample_id"] = config["sample_experiment"][wildcards.fastq]
    readgroup = '-r ID:{s_id}_{f_id}_{f_lane} -r SM:{s_id} -r LB:{lib} -r PU:{f_lane} -r PL:ILLUMINA'.format(f_id=samples["flowcell_id"],
                                                                                                             f_lane=samples["flowcell_lane"],
                                                                                                             lib=samples["library"],
                                                                                                             s_id=samples["sample_id"])
    #print(readgroup)
    return readgroup

rule add_readgroup:
    input: 
         'output/Omni-C_pairsam/{species}/{ccs_settings}/{hifi_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam'
    output: 
         'output/Omni-C_pairsam/{species}/{ccs_settings}/{hifi_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.rg.bam'
    params: 
        readgroup = get_seqinfo
    conda: "../envs/STAR-EBSeq-RSEM.yaml"
    shell:
        'samtools addreplacerg {params.readgroup} -o {output} {input}'


rule index_bams_rg:
    input:
         'output/Omni-C_pairsam/{species}/{ccs_settings}/{hifi_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.rg.bam'
    output:
         'output/Omni-C_pairsam/{species}/{ccs_settings}/{hifi_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.rg.bam.bai'
    threads: 40
    params:
        slurm_opts=lambda wildcards: "-n1 "
                                     "--share "
                                     "--export ALL "
                                     "--mem 30000 "
                                     "--time 0-6:00:00 "
                                     "-J index_{sample} "
                                     "-o logs/index_{sample}_%j.logs "
                                     "-p defq "
                                     "".format(sample=wildcards.fastq)
    shell:
        "samtools index -@ {threads} {input}"


rule pretext_symlink_genome:
    input: 
        fasta = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa"
    output: 
        fasta = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/ref.fa",
    conda: "../envs/RapidCuration.yaml"
    threads: 2
    shell: "ROOTDIR=$(pwd -P); "
           "ln -sf $ROOTDIR/{input.fasta} {output.fasta}; "


rule pretext_symlinks_bam:
    input: 
        bam = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam",
    output: 
        #bam = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.rg.bam",
        bam = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/bam.mapq{mapq}.sorted.rg.bam",
    conda: "../envs/RapidCuration.yaml"
    threads: 2
    shell: "ROOTDIR=$(pwd -P); "
           "ln -sf $ROOTDIR/{input.bam} {output.bam}; "


rule pretext_symlinks_fofn:
    input: 
        bam = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam",
    output: 
        bam = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/fofn.mapq{mapq}.cram",
    params:
        outdir = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/"
    conda: "../envs/RapidCuration.yaml"
    threads: 2
    shell: "ROOTDIR=$(pwd -P); "
           "echo '{wildcards.species}.{wildcards.genometype}.{wildcards.hic}.mapq{wildcards.mapq}.sorted.bam' > {output}"
           "ln -sf fofn.mapq{mapq}.cram {params.outdir}/fofn.cram"


rule pretext_map:
    input: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/bam.mapq{mapq}.sorted.rg.bam"
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pretext.mapq{mapq}.pretext"
    conda: "../envs/RapidCuration.yaml"
    #params:
    #    outdir = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output"
    threads: 2
    shell: "samtools view -h {input} | "
           " awk '{{print gensub(/scaff_([0123456789]+)_contig_tig([0123456789]+)[^\t]+(\s)/,\"\\1_\\2\\3\",\"g\")}}' | "
           " PretextMap -o {output} --sortby 'length' --sortorder 'descend' --highRes"

rule rapidcuration_Gap:
    input: 
        #bam = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam",
        fasta = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/ref.fa",
        #fofn = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/cram.fofn",
    #output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}.mapq{mapq}_gap.bedgraph",
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/gap.mapq{mapq}.done",
    params:
        sample = "{species}.{genometype}.{hic}",
        base_dir1 = "/global/scratch", 
        base_dir2 = "/global/scratch2", 
        base_dir3 = "/global/home/users/mvazquez",
        nfs_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly",
        lustre_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly",
        data_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly/output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/",
        output_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly/output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output"
    conda: "../envs/RapidCuration.yaml"
    shell: "ROOTDIR=$(pwd -P); "
           "cd {params.output_dir} && "
           "singularity run "
           "-B {params.base_dir1} "
           "-B {params.base_dir2} "
           "-B {params.base_dir3} "
           "-B {params.nfs_dir}:/nfs "
           "-B {params.lustre_dir}:/lustre "
           "-B {params.data_dir}:/data "
           "-B {params.output_dir}:/output "
           "$ROOTDIR/code/rapid-curation/rapid_hic_software/runGap.sif -t {params.sample}; "
           "touch {output}"

rule pt_addGap:
    input: 
        gap = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}.mapq{mapq}_gap.bedgraph",
        pt = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}.mapq{mapq}.pretext"
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pt_addGap.done"
    conda: "../envs/RapidCuration.yaml"
    shell: "cat {input.gap} | PretextGraph -i {input.pt} -n 'Gap' && touch {output}"


rule rapidcuration_telo:
    input: 
        #bam = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam",
        fasta = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/ref.fa",
        #fofn = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/cram.fofn",
    #output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}.mapq{mapq}_gap.bedgraph_telomere.bedgraph",
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/telo.done",
    params:
        sample = "{species}.{genometype}.{hic}",
        base_dir1 = "/global/scratch", 
        base_dir2 = "/global/scratch2", 
        base_dir3 = "/global/home/users/mvazquez",
        nfs_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly",
        lustre_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly",
        data_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly/output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}",
        output_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly/output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output"
    conda: "../envs/RapidCuration.yaml"
    shell: "ROOTDIR=$(pwd -P); "
           "cd {params.output_dir} && "
           "singularity run "
           "-B {params.base_dir1} "
           "-B {params.base_dir2} "
           "-B {params.base_dir3} "
           "-B {params.nfs_dir}:/nfs "
           "-B {params.lustre_dir}:/lustre "
           "-B {params.data_dir}:/data "
           "-B {params.output_dir}:/output "
           "$ROOTDIR/code/rapid-curation/rapid_hic_software/runTelo.sif -s {params.sample} -t {params.sample}; "
           "touch {output}"


rule pt_addTelo:
    input: 
        telo = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}_gap.bedgraph_telomere.bedgraph",
        pt = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}.mapq{mapq}.pretext"
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pt_addTelo.mapq{mapq}.done"
    conda: "../envs/RapidCuration.yaml"
    shell: "cat {input.telo} | PretextGraph -i {input.pt} -n 'Telomere' && touch {output}"


rule rapidcuration_Repeat:
    input: 
        #bam = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/{species}.{genometype}.{hic}.mapq{mapq}.sorted.bam",
        fasta = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/ref.fa",
        #fofn = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/cram.fofn",
    #output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}_repeat_density.mapq{mapq}.bw",
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/repeat.mapq{mapq}.done",
    params:
        sample = "{species}.{genometype}.{hic}",
        base_dir1 = "/global/scratch", 
        base_dir2 = "/global/scratch2", 
        base_dir3 = "/global/home/users/mvazquez",
        nfs_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly",
        lustre_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly",
        data_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly/output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}",
        output_dir = "/global/scratch2/mvazquez/projects/BatGenomeAssembly/output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output"
    conda: "../envs/RapidCuration.yaml"
    shell: "ROOTDIR=$(pwd -P); "
           "cd {params.output_dir} && "
           "singularity run "
           "-B {params.base_dir1} "
           "-B {params.base_dir2} "
           "-B {params.base_dir3} "
           "-B {params.nfs_dir}:/nfs "
           "-B {params.lustre_dir}:/lustre "
           "-B {params.data_dir}:/data "
           "-B {params.output_dir}:/output "
           "$ROOTDIR/code/rapid-curation/rapid_hic_software/runRepeat.sif -t {params.sample}; "
           "touch {output}"


rule pt_addRep:
    input: 
        rep_bw = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}.mapq{mapq}_repeat_density.bw",
        pt = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/{species}.{genometype}.{hic}.mapq{mapq}.pretext"
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pt_addRep.mapq{mapq}.done"
    conda: "../envs/RapidCuration.yaml"
    shell: "bigWigToBedGraph {input.rep_bw} /dev/stdout | PretextGraph -i {input.pt} -n 'Repeats' && touch {output}"


rule pt_alladded:
    input: ["output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pt_addRep.done", "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pt_addTelo.done", "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pt_addGap.done"]
    output: "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/pt_addall.done"
    shell: "touch {output}"



rule pt_rapid_join:
    input: 
        fa = "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa",
        pt_tpf = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/rapid_prtxt_XL.tpf",
        pt_chrs = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output/chrs.csv",
    output: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}_manualCuration.fa"
    conda: "../envs/RapidCuration.yaml"
    params:
        pt_dir = "output/pretext_map/{species}/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}/output"
    shell: "PROJDIR=$(pwd -P); "
           "cd {params.pt_dir}; "
           "/global/home/users/mvazquez/miniconda3/envs/RapidCuration/bin/perl $PROJDIR/code/rapid-curation/rapid_join.pl -fa $PROJDIR/{input.fa} -tpf $PROJDIR/{input.pt_tpf} -csv $PROJDIR/{input.pt_chrs} -out $PROJDIR/{output}; "
           "cd - "
