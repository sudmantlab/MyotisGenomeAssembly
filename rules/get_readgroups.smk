configfile: "config_trimmed.json"

def get_all_fastq(wildcards):
    inputs = []
    samples = list(config['sample_experiment'].keys())
    for sample in samples:
        for fastq in config['input_samples'][sample]:
            inputs.append("{sample_path}/"
                           "{sample_fastq}"
                           "".format(sample_path=config['input_sample_path'],
                                     sample_fastq=fastq))
    return inputs

rule get_readgroupinfo:
    input: get_all_fastq
    output: 'data/RNA-seq/seqinfo.tsv'
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
