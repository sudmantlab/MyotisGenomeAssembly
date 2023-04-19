rule funannotate_docker_clean:
    input: "data/genomes/{genome}.fa"
    output: "data/genomes/{genome}.cleaned.fa"
    shell: "./funannotate-docker clean -i $(pwd -P)/{input} -o $(pwd -P)/{output} -m 10000"
