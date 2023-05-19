rule trimmomatic_pe:
    input:
        r1="data/RNA-seq/{species}/{sample}_R1_001.fastq.gz",
        r2="data/RNA-seq/{species}/{sample}_R2_001.fastq.gz"
    output:
        r1="output/trimmed/{species}/{sample}-trimmed_R1.fastq.gz",
        r2="output/trimmed/{species}/{sample}-trimmed_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{species}/{sample}-trimmed_R1_unpaired.fastq.gz",
        r2_unpaired="output/trimmed/{species}/{sample}-trimmed_R2_unpaired.fastq.gz"
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
    wrapper:
        "0.74.0/bio/trimmomatic/pe"


rule trimmomatic_pe_10x:
    input:
        r1="data/10x/{species}/{sample}_R1_001.fastq.gz",
        r2="data/10x/{species}/{sample}_R2_001.fastq.gz"
    output:
        r1="output/trimmed-10x/{species}/{sample}-trimmed_R1.fastq.gz",
        r2="output/trimmed-10x/{species}/{sample}-trimmed_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{species}/{sample}-trimmed_R1_unpaired.fastq.gz",
        r2_unpaired="output/trimmed/{species}/{sample}-trimmed_R2_unpaired.fastq.gz"
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
    threads: 20
    wrapper:
        "0.74.0/bio/trimmomatic/pe"


rule trimmomatic_pe_HiC:
    input:
        r1="data/Hi-C/{species}/{sample}_R1_001.fastq.gz",
        r2="data/Hi-C/{species}/{sample}_R2_001.fastq.gz"
    output:
        r1="output/trimmed-hic/{species}/{sample}-trimmed_R1.fastq.gz",
        r2="output/trimmed-hic/{species}/{sample}-trimmed_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{species}/{sample}-trimmed_R1_unpaired.fastq.gz",
        r2_unpaired="output/trimmed/{species}/{sample}-trimmed_R2_unpaired.fastq.gz"
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
    threads: 32
    wrapper:
        "0.74.0/bio/trimmomatic/pe"
