
def get_ccs_stats(wildcards):
    hifi_path = "output/HiFi-adapterFiltered/{species}/{settings}/{pacbio1}/{pacbio2}/{id}.ccs.filt.fastq.gz"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Species"] == wildcards.species]
    samples = samples.to_records(index=False)
    input_samples = [hifi_path.format(species=s[0], settings = wildcards.settings, pacbio1 = s[1], pacbio2 = s[2], id = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for species {}. Check pepsamples.tsv and try again!".format(wildcards.species))
    else:
        return input_samples

#rule all_readsDepth_CCS:
#    version: "0.1"
#    input: lambda wildcards: 


rule get_readsDepth_CCS:
    version: "0.1"
    input: get_hifiasm_inputs
    output: "output/readDepthStats_CCS/{settings}/{species}_CCS_stats.tsv"
    run:
        import json
        from pathlib import Path
        outlist = []
        for report in input:
            with Path(report).open() as f:
                data = json.load(f)
            read = [i['value'] for i in data["attributes"] if i['id'] == "ccs2.number_of_ccs_reads" ]
            bases = [i['value'] for i in data["attributes"] if i['id'] == "ccs2.total_number_of_ccs_bases" ]
            name = report.stem
            outlist.append("\t".join([name, str(read[0]), str(bases[0])]) + "\n")
        with Path(output).open("w") as f:
            f.writelines(outlist)


rule get_hifiasm_merqury_stats:

