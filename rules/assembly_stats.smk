#rule run_assembly_stats:
#    input: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.fa"
#    output: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.assembly_stats"
#    shell: "code/assemblystats/target/release/assemblystats {input} > {output}"


#rule get_Nx:
#    input: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.assembly_stats"
#    output: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.Nx.csv"
#    shell: "echo Nx,{wildcards.species}.{wildcards.genometype}.{wildcards.hic} > {output}; "
#           "cat {input} | grep -oe 'N[0-9]\+: [0-9]*' | sed 's/: /,/' >> {output}"

#rule get_Lx:
#    input: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.assembly_stats"
#    output: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{genometype}.{hic}.Lx.csv"
#    shell: "echo Lx,{wildcards.species}.{wildcards.genometype}.{wildcards.hic} > {output}; "
#           "cat {input} | grep -oe 'L[0-9]\+: [0-9]*' | sed 's/: /,/' >> {output}"


rule run_assembly_stats:
    input: "data/assemblies/{ccs_opts}/{hifiasm_opts}/{species}.{assembly}.fa"
    output: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{assembly}.assembly_stats"
    shell: "code/assemblystats/target/release/assemblystats {input} > {output}"


rule get_Nx:
    input: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{assembly}.assembly_stats"
    output: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{assembly}.Nx.csv"
    shell: "echo Nx,{wildcards.species}.{wildcards.assembly} > {output}; "
           "cat {input} | grep -oe 'N[0-9]\+: [0-9]*' | sed 's/: /,/' >> {output}"

rule get_Lx:
    input: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{assembly}.assembly_stats"
    output: "output/assembly_stats/{ccs_opts}/{hifiasm_opts}/{species}.{assembly}.Lx.csv"
    shell: "echo Lx,{wildcards.species}.{wildcards.assembly} > {output}; "
           "cat {input} | grep -oe 'L[0-9]\+: [0-9]*' | sed 's/: /,/' >> {output}"

   
