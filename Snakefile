# Genome Assembly
include: "rules/HiFi.smk"

# Genome QC
include: "rules/genomeQC.smk"

# Trinity
include: "rules/QC_RNASeq.smk"
include: "rules/trimmomatic.smk"
include: "rules/trinity.smk"

# Annotate
include: "rules/mitoHiFi.smk"

# Omni-C mapping to genome
include: "rules/omniC.smk"
