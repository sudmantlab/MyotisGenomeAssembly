rule tidk_search:
    input: "data/genomes/{genome}.fa"
    output: multiext("output/tidk/{genome}/{genome}_telo_telomeric_repeat_windows", ".bedgraph", ".csv")
    threads: 1
    shell: "tidk search --fasta {input} -d output/tidk/{wildcards.genome} -o {wildcards.genome}_telo -e bedgraph --string AACCCT; "
           "tidk search --fasta {input} -d output/tidk/{wildcards.genome} -o {wildcards.genome}_telo -e csv --string AACCCT"


rule tidk_plot:
    input: "output/tidk/{genome}/{genome}_telo_telomeric_repeat_windows.csv"
    output: "output/tidk/plots/{genome}_telo.svg"
    threads: 1
    shell: "tidk plot -c {input} -o output/tidk/plots/{wildcards.genome}_telo"

