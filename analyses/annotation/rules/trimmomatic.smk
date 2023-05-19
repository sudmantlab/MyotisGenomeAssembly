rule trimmomatic_pe:
    input:
        r1="data/RNA-seq/{species}/PAIRED/{sample}_R1.fastq.gz",
        r2="data/RNA-seq/{species}/PAIRED/{sample}_R2.fastq.gz"
    output:
        r1="output/trimmed/{species}/PAIRED/{sample}-trimmed_R1.fastq.gz",
        r2="output/trimmed/{species}/PAIRED/{sample}-trimmed_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{species}/PAIRED/{sample}-trimmed_R1_unpaired.fastq.gz",
        r2_unpaired="output/trimmed/{species}/PAIRED/{sample}-trimmed_R2_unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{species}/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:data/trimmomatic-adapters/TruSeq3-PE-2.fa:2:40:15", 
                 "SLIDINGWINDOW:5:20"
                ],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: 4
    resources:
        tmpdir = "/global/scratch2/mvazquez/tmp"
    wrapper:
        "0.74.0/bio/trimmomatic/pe"

rule trimmomatic_se:
    input:
        "data/RNA-seq/{species}/SINGLE/{sample}.fastq.gz"
    output:
        "output/trimmed/{species}/SINGLE/{sample}-trimmed.fastq.gz",
    log:
        "logs/trimmomatic/{species}/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:data/trimmomatic-adapters/TruSeq3-PE-2.fa:2:40:15", 
                 "SLIDINGWINDOW:5:20"
                ],
        # optional parameters
        extra="",
        compression_level="-9"
    resources:
        tmpdir = "/global/scratch2/mvazquez/tmp"
    threads: 4
    wrapper:
        "0.74.0/bio/trimmomatic/se"
