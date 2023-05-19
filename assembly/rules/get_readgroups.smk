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
