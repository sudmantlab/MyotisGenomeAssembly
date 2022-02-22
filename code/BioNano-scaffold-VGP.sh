#!/bin/bash

## Modified version of https://raw.githubusercontent.com/VGP/vgp-assembly/master/pipeline/bionano/hybrid_scaffold.sh

name=$1
ENZYME=$2
BMAP=../data/BMAP/$ENZYME.cmap # DLE1.cmap
ASM=$3	# ln -s to the asm.fasta
CONFIG=$4
OUTPUT=$5
RefAligner=$tools/bionano/Solve3.6.1/RefAligner/7437.7523rel/avx/RefAligner

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "\
perl $tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/hybridScaffold.pl \
        -n $ASM \
	-b $BMAP \
	-c $CONFIG \
	-r $RefAligner \
	-B 2 \
	-N 2 \
	-f \
        -o $OUTPUT
"
perl ../code/bionano/Solve3.6.1/HybridScaffold/04122018/hybridScaffold.pl \
        -n $ASM \
	-b $BMAP \
	-c $CONFIG \
	-r $RefAligner \
	-B 2 \
	-N 2 \
	-f \
        -o $OUTPUT
