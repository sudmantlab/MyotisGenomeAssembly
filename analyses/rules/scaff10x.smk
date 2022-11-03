rule scaff10x_hifi:
    input: "hifi/genome/path"
    output: "10x/genome/path"
    threads: 30
    shell: """
src/scaff10x \
-nodes {threads} \
-longread 1 \
-gap 100 \
-matrix 2000 \
-reads 10 \
-score 20 \
-edge 50000 \
-link 8 \
-block 50000 \
-plot barcode_lengtg.png \
<input_assembly_fasta/q_file> \
<Input_read_1> \ 
<Input_read_2> \
<Output_scaffold_file>
"""
