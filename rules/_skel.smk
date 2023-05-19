rule run_assembly_stats:
    version: 0.1
    input: "data/genomes/{genome}.fa"
    output: temp("output/assembly_stats/{genome}.assembly_stats")
    shell: "code/assemblystats/target/release/assemblystats {input} > {output}"

   
