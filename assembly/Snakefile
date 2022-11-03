configfile: "config.yaml"

# Genome Assembly
include: "rules/HiFi.smk"

# Genome QC
include: "rules/genomeQC.smk"
include: "rules/summaryStatTables.smk"


# Trinity
include: "rules/QC_RNASeq.smk"
include: "rules/trimmomatic.smk"
include: "rules/trinity.smk"

# Annotate
include: "rules/mitoHiFi.smk"

# Omni-C mapping to genome
include: "rules/omniC.smk"
include: "rules/omniC-scaffold.smk"

# Scaffolding
include: "rules/yahs.smk"

# Curation
include: "rules/cooler.smk"
include: "rules/pretext.smk"


# Experimental
include: "rules/deepConsensus.smk"
include: "rules/quast.smk"
