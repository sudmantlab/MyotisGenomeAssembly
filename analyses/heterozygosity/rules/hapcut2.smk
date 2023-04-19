rule extractHAIRS: 
	input: 
		bam = "output/minimap2-split/{species}/{contig}.sorted.bam",
		vcf = "output/freebayes-calls-split/M_lucifugus-ptg_4.vcf"
	output: "output/hapcut2-HAIRS/{species}/{contig}.fragments"
	conda: "../envs/hapcut2.yml"
	shell: "extractHAIRS --pacbio 1 --indels --bam {input.bam} --VCF {input.vcf} --out {output}"

rule 
	input: 
		vcf = "output/freebayes-calls-split/M_lucifugus-ptg_4.vcf",
		fragments = "output/hapcut2-HAIRS/{species}/{contig}.fragments"
	output: "output/hapcut2/{species}/{species}-{contig}.phased.block"
	conda: "../envs/hapcut2.yml"
	shell:" HAPCUT2 --ea 1 --fragments {input.fragments} --VCF {input.vcf} --output haplotype_output_file"
