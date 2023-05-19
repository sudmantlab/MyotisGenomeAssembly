configfile: "config.yaml"

localrules: link_hifiasm_assembly, link_yahs_assembly, pretext_symlinks, symlink_assemblies, add_readgroup, index_bams_rg

wildcard_constraints:
    # binSize_bp = get_binsize_bp,  ## rules/cooler.smk
    # chunk="\d{1,2}",  ## rules/CCS_split.smk
    isDefault="(defaults-)?",
    genometype = "(hap[12].)?[pr]_[cu]tg",
    contigtype = "[pr]_[cu]tg",
    species="[A-Z]_[a-z]+",
    hic = "(bp)|(hic)(-purgedDups)?(_mapq[0-9]+)?(_noContigEC_minQ[0-9]+)?(_scaffolds_final)?(-purgedDups)?(_mapq[0-9]+)?(_manualCuration)?",
    mapq = "[0-9]+",
    minmapq = "[0-9]+",
    unit = "bp|kb|mb"


# Genome Assembly
include: "rules/HiFi.smk"
include: "rules/purge_dups.smk"

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
#include: "rules/omniC-scaffold.smk"

# Scaffolding
include: "rules/yahs.smk"

# Curation
include: "rules/cooler.smk"
include: "rules/pretext.smk"
include: "rules/genomes_curate.smk"


# Experimental
include: "rules/deepConsensus.smk"
include: "rules/quast.smk"
include: "rules/ragtag.smk"

