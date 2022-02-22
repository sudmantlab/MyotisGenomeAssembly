rule bionano_scaffold_vgp_fromHiFi:
    input:
        BMAP = "data/BMAP/{enzyme}.cmap",
        assemblyFasta = "output/hifiasm-fasta/{species}/{settings}/{genome_opts}/{species}.p_ctg.fa",
        configXML = "data/bionano_config/bionanoConfig.xml"
    output: directory("output/bioNano_scaffold/{species}/bioNano_{enzyme}/{settings}/{genome_opts}/{species}.p_ctg.bnsc.fa")
    shell: "../code/BioNano-scaffold-VGP.sh {wildcards.genomeName} {wildcards.enzyme} {input.assemblyFasta} {input.configXML}"


rule bionano_scaffold_vgp_from10x:
    input:
        BMAP = "data/BMAP/{enzyme}.cmap",
        assemblyFasta = "output/hifiasm-fasta/{species}/{settings}/{genome_opts}/{species}.p_ctg.fa",
        configXML = "data/bionano_config/bionanoConfig.xml"
    output: directory("output/bioNano_scaffold/{species}/bioNano_{enzyme}/{settings}/{genome_opts}/{species}.p_ctg.bnsc.fa")
    shell: "../code/BioNano-scaffold-VGP.sh {wildcards.genomeName} {wildcards.enzyme} {input.assemblyFasta} {input.configXML}"
	
