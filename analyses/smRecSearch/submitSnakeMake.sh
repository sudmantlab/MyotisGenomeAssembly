#!/bin/bash

conda activate R-GeneLists

snakemake \
    -kpr \
    --ri \
    --nolock \
    -j 500 \
    --cluster-config /project2/gilad/juanvazquez/projects/smRecSearch/cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --ntasks-per-node={cluster.tasks} \
        --partition=broadwl \
        " \
    $*


#    --restart-times 2 \
#    --ri \

