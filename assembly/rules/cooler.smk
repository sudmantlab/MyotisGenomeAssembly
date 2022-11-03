rule bgzip_pairs:
    input: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs"
    output: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs.gz"
    conda: "../envs/RapidCuration.yaml"
    shell: "bgzip {input}"


rule pairix:
    input: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs.gz"
    output: "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs.gz.px2"
    conda: "../envs/RapidCuration.yaml"
    shell: "pairix {input}"


def get_binsize_bp(wildcards):
    units = wildcards.unit
    size = wildcards.binSize
    if units.lower() == "bp":
        multiplier = 1
    elif units.lower() == "kb":
        multiplier = 1000
    elif units.lower() == "mb":
        multiplier = 1000000
    else:
        raise Exception("Not a real unit!")
    return str(int(size) * multiplier)


rule cooler_matrix:
    input: 
        pairs = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs.gz",
        idx = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs.gz.px2",
        genome = "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.genome"
    output: "output/cooler/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.matrix_{binSize}-{unit}.cool"
    wildcard_constraints:
        unit = "bp|kb|mb"
    params:
        binSize_bp = get_binsize_bp
    conda: "../envs/RapidCuration.yaml"
    threads: 52
    shell: "cooler cload pairix -p {threads} {input.genome}:{params.binSize_bp} {input.pairs} {output}"


rule cooler_matrix_yahs:
    input: 
        pairs = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs.gz",
        idx = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.pairs.gz.px2",
        genome = "output/yahs/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.genome"
    output: "output/cooler/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.matrix_{binSize}-{unit}.cool"
    wildcard_constraints:
        unit = "bp|kb|mb"
    params:
        binSize_bp = get_binsize_bp
    conda: "../envs/RapidCuration.yaml"
    threads: 52
    shell: "cooler cload pairix -p {threads} {input.genome}:{params.binSize_bp} {input.pairs} {output}"


rule cooler_zoomify:
    input: "output/cooler/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.matrix_{binSize}-{unit}.cool"
    output: "output/cooler/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.{hic}.matrix_{binSize}-{unit}.mcool"
    conda: "../envs/RapidCuration.yaml"
    threads: 52
    shell: "cooler zoomify -r N --balance -p {threads} {input}"
