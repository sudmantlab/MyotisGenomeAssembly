#rule pipe_zcat:
#    input: "data/genomes/{genome}.fa.gz"
#    output: pipe("data/genomes/{genome}.fa")
#    shell: "zcat {input} > {output}"

rule run_assembly_stats:
    version: 0.2
    input: "data/genomes/{genome}.fa"
    output: temp("output/assembly_stats/{genome}.assembly_stats")
    shell: "code/assemblystats/target/release/assemblystats --genomename {wildcards.genome} {input} {output}"


#rule get_Nx:
#    version: 0.1
#    input: "output/assembly_stats/{genome}.assembly_stats"
#    output: "output/assembly_stats/{genome}.Nx.csv"
#    shell: "echo Nx,{wildcards.genome} > {output}; "
#           "cat {input} | grep -oe 'N[0-9]\+: [0-9]*' | sed 's/: /,/' >> {output}"

#rule get_Lx:
#    version: 0.1
#    input: "output/assembly_stats/{genome}.assembly_stats"
#    output: "output/assembly_stats/{genome}.Lx.csv"
#    shell: "echo Lx,{wildcards.genome} > {output}; "
#           "cat {input} | grep -oe 'L[0-9]\+: [0-9]*' | sed 's/: /,/' >> {output}"

rule all_summary_stats:
    input: expand("output/assembly_stats/{genome}.assembly_stats", genome=glob_wildcards('data/genomes/{genome}.fa')[0])
    output: "output/assembly_stats/chiroptera.assembly_stats"
    shell:
        "echo -e 'genomeName\tstat\tv1\v2\n' > {output}; "
        "cat {input} | "
        " grep -v 'contigs >=' | "
        " grep -v 'genomeName' >> {output}"
        
