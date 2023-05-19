rule diamond_makedb:
    input: 
    output: 
    params:
        taxonnodes = "/global/scratch/users/mvazquez/taxdump/nodes.dmp",
        taxonnames = "/global/scratch/users/mvazquez/taxdump/names.dmp"
    #shell: "diamond makedb --in {input} -d {wildcards.genome} --taxonnodes {params.taxonnodes} --taxonnames {params.taxonnames}"
    shell: "PROJDIR=$(pwd -P);"
           "cd data/
           "diamond makedb --in {input} -d {wildcards.genome}"
