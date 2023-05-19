rule trim_reads_PE:
  version: "0.9" # will become 1.0 once I move params to output path instead of hardcoding them.
  input: 
    forward = "data/RNASeq/{run}/{sample}_R1_001.fastq.gz",
    reverse = "data/RNASeq/{run}/{sample}_R2_001.fastq.gz"
  output: 
    forward_paired = "output/trimmomatic/{run}/{sample}_R1_001_paired.fastq.gz",
    forward_unpaired = "output/trimmomatic/{run}/{sample}_R1_001_unpaired.fastq.gz",
    reverse_paired = "output/trimmomatic/{run}/{sample}_R2_001_paired.fastq.gz",
    reverse_unpaired = "output/trimmomatic/{run}/{sample}_R2_001_unpaired.fastq.gz"
  params:
    illuminaclip = "TruSeq3-PE.fa:2:30:10:2:keepBothReads",
    # leading = "3",
    # trailing = "3",
    slidingwindow = "5:20",
    minlen = "75"
 threads: 32
 conda: "envs/readTrimming.yml"
 shell:
   "trimmomatic PE -threads {threads}"
   "{input.forward} {input.reverse} "
   "{output.forward_paired} {output.forward_unpaired} "
   "{output.reverse_paired} {output.reverse_unpaired} "
   " ILLUMINACLIP:{params.illuminaclip} "
   # "LEADING:{params.leading} "
   # "TRAILING:{params.trailing} "
   "SLIDINGWINDOW:{params.slidingwindow} "
   "MINLEN:{params.minlen} "
